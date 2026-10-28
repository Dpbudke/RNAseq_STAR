#!/usr/bin/env python3

import os
import sys
import gzip
from collections import Counter
from Bio import SeqIO
from scipy.stats import spearmanr
import csv

def get_kmer_profile(fastq_file, k=7, n_reads=100000):
    """
    Generate a k-mer frequency profile from the first n reads of a FASTQ file.
    
    Args:
        fastq_file (str): Path to FASTQ file (can be gzipped)
        k (int): k-mer size
        n_reads (int): Number of reads to process
    
    Returns:
        dict: k-mer frequencies
    """
    kmer_counts = Counter()
    
    # Handle both gzipped and regular files
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'
    
    try:
        with open_func(fastq_file, mode) as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                if i >= n_reads:
                    break
                # Generate k-mers from sequence
                seq = str(record.seq)
                for j in range(len(seq) - k + 1):
                    kmer = seq[j:j+k]
                    kmer_counts[kmer] += 1
    except FileNotFoundError:
        print(f"Warning: Could not find file {fastq_file}")
        return None
    except Exception as e:
        print(f"Error processing file {fastq_file}: {str(e)}")
        return None
    
    # Normalize frequencies
    total = sum(kmer_counts.values())
    if total == 0:
        return None
    kmer_freqs = {k: v/total for k, v in kmer_counts.items()}
    return kmer_freqs

def compare_samples(raw_dir, proc_dir, samples, k=7, layout="paired"):
    """
    Compare k-mer profiles between raw and processed samples for both R1 and R2.
    
    Args:
        raw_dir (str): Path to raw reads directory
        proc_dir (str): Path to processed reads directory
        samples (list): List of sample names
        k (int): k-mer size for profiling
    
    Returns:
        tuple: Dictionaries containing correlation values for R1 and R2
    """
    # Store k-mer profiles
    profiles_r1 = {}
    profiles_r2 = {}
    
    for sample in samples:
        print(f"\nProcessing sample {sample}...")
        
        # Define file paths
        raw_r1 = os.path.join(raw_dir, sample, f"{sample}_R1_concatenated.fastq.gz")
        raw_r2 = os.path.join(raw_dir, sample, f"{sample}_R2_concatenated.fastq.gz")
        # Single-end runs: the trimmed reads are <sample>_SE.fastq.gz and there is no R2.
        proc_r1 = os.path.join(proc_dir, sample, f"{sample}_R1.fastq.gz" if layout == "paired" else f"{sample}_SE.fastq.gz")
        proc_r2 = os.path.join(proc_dir, sample, f"{sample}_R2.fastq.gz")
        
        # Generate profiles for R1
        print(f"Processing R1 files...")
        raw_profile_r1 = get_kmer_profile(raw_r1)
        proc_profile_r1 = get_kmer_profile(proc_r1)
        
        if raw_profile_r1 and proc_profile_r1:
            profiles_r1[f"{sample}_raw"] = raw_profile_r1
            profiles_r1[f"{sample}_proc"] = proc_profile_r1
        
        if layout != "paired":
            continue

        # Generate profiles for R2
        print(f"Processing R2 files...")
        raw_profile_r2 = get_kmer_profile(raw_r2)
        proc_profile_r2 = get_kmer_profile(proc_r2)
        
        if raw_profile_r2 and proc_profile_r2:
            profiles_r2[f"{sample}_raw"] = raw_profile_r2
            profiles_r2[f"{sample}_proc"] = proc_profile_r2
    
    return profiles_r1, profiles_r2

def calculate_correlations(profiles):
    """
    Calculate correlation matrix for a set of profiles.
    
    Args:
        profiles (dict): Dictionary of k-mer profiles
    
    Returns:
        dict: Dictionary of correlation values
    """
    correlations = {}
    for s1 in profiles.keys():
        correlations[s1] = {}
        kmers = sorted(set(profiles[s1].keys()))
        s1_values = [profiles[s1].get(k, 0) for k in kmers]
        
        for s2 in profiles.keys():
            s2_values = [profiles[s2].get(k, 0) for k in kmers]
            corr, _ = spearmanr(s1_values, s2_values)
            correlations[s1][s2] = corr
    
    return correlations

def write_correlations(correlations, output_file):
    """
    Write correlation results to a CSV file.
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        samples = sorted(correlations.keys())
        writer.writerow(['Sample'] + samples)
        for s1 in samples:
            row = [s1] + [correlations[s1][s2] for s2 in samples]
            writer.writerow(row)

def print_correlation_summary(correlations, samples, read_type):
    """
    Print correlation summary including matched and mismatched samples.
    """
    print(f"\n{read_type} Sample Correlation Summary:")
    print("\nMatched Sample Correlations (should be high):")
    for sample in samples:
        raw_key = f"{sample}_raw"
        proc_key = f"{sample}_proc"
        if raw_key in correlations and proc_key in correlations[raw_key]:
            corr = correlations[raw_key][proc_key]
            print(f"{sample}: Raw vs Processed correlation = {corr:.3f}")
    
    print("\nMismatch Sample Correlations (should be lower):")
    for i, sample1 in enumerate(samples):
        for sample2 in samples[i+1:]:
            raw_key = f"{sample1}_raw"
            proc_key = f"{sample2}_proc"
            if raw_key in correlations and proc_key in correlations[raw_key]:
                corr = correlations[raw_key][proc_key]
                print(f"{sample1} raw vs {sample2} processed correlation = {corr:.3f}")

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="k-mer profile correlation between raw and HTStream-processed reads")
    ap.add_argument("--raw-dir", default=os.environ.get("CONCAT_DIR", "concatenated_reads"))
    ap.add_argument("--proc-dir", default=os.environ.get("HTS_OUT", "01-HTS_Preproc"))
    ap.add_argument("--samples", default=os.environ.get("SAMPLES_TXT", "samples.txt"),
                    help="file with one sample per line")
    ap.add_argument("--n", type=int, default=5,
                    help="how many samples to take from the top of that file (0 = all)")
    ap.add_argument("--out-prefix", default="sample_correlations")
    ap.add_argument("--layout", choices=["paired", "single"], default=os.environ.get("LIBRARY_LAYOUT", "paired"))
    a = ap.parse_args()

    if not os.path.isfile(a.samples):
        sys.exit(f"sample list not found: {a.samples}")
    with open(a.samples) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    if a.n > 0:
        samples = samples[:a.n]
    if not samples:
        sys.exit(f"{a.samples} is empty")

    print(f"raw:  {a.raw_dir}")
    print(f"proc: {a.proc_dir}")
    print(f"samples: {', '.join(samples)}")

    profiles_r1, profiles_r2 = compare_samples(a.raw_dir, a.proc_dir, samples, layout=a.layout)
    if not profiles_r1 and not profiles_r2:
        sys.exit("no profiles were built. check --raw-dir and --proc-dir.")

    correlations_r1 = calculate_correlations(profiles_r1)
    write_correlations(correlations_r1, f"{a.out_prefix}_R1.csv")
    if a.layout == "paired":
        correlations_r2 = calculate_correlations(profiles_r2)
        write_correlations(correlations_r2, f"{a.out_prefix}_R2.csv")
        print(f"\nwrote {a.out_prefix}_R1.csv and {a.out_prefix}_R2.csv")
    else:
        print(f"\nwrote {a.out_prefix}_R1.csv (single-end)")

    print_correlation_summary(correlations_r1, samples, "R1")
    if a.layout == "paired":
        print_correlation_summary(correlations_r2, samples, "R2")


if __name__ == "__main__":
    main()
