# Changelog

## 2026-10-28

- References for other organisms: `REFERENCE_SOURCE` selects GENCODE mouse (default, unchanged), GENCODE human,
  any Ensembl species and release, or your own FASTA and GTF. Downloads are checked against the provider's
  checksums, SHA-256 is recorded, and the genome and annotation must share chromosome names. For other sources the
  HTStream rRNA count uses the annotation's rRNA genes or `RRNA_FASTA`. Small genomes get a smaller STAR
  suffix-array index; mouse and human indexes are unchanged.
- Single-end libraries (`LIBRARY_LAYOUT=single`), lanes without a provider checksum file (`SOURCE_MD5=none`
  records SHA-256 instead) and a configurable post-trim minimum length (`HTS_MIN_LENGTH`, default 50; `check`
  refuses reads shorter than it).
- `SOFTWARE_MODE=container`: STAR, SAMtools, HTStream and the QC Python libraries from a pinned image
  (`container/Dockerfile`, built by `.github/workflows/container.yml` as `ghcr.io/bjbennettlab/rnaseq_star`). On
  a Bennett lab DO sample (paired-end) and the yeast single-end set it gave the same results as the Ceres modules; it carries the module's exact HTStream build.
- `HTS_DEDUP=no` skips HTStream duplicate removal (default `yes`, unchanged). On the yeast single-end test it
  removed 83 to 85% of reads; `check` now warns when it is on for single-end reads.
- `check` refuses a Slurm account the user does not have and work or reference folders outside that SCINet
  project's `/project` or `/90daydata`.
- Checked on Ceres (October 28): a Bennett lab DO sample (two lanes, paired-end, GENCODE mouse) through this version
  gave the same gene counts, mapping statistics, trimmed reads, HTStream counts and reference manifests as the
  September 10, 2026 run; a public yeast single-end set (ENA PRJNA1365179, Ensembl release 113) ran through every
  stage and passed validation.

## 2026-10-26

- With `EXECUTOR=slurm`, the stages usually run by hand (`stats`, `counts`, `htsqc`, `confirm`, `kmerqc`,
  `compare`) now submit a small Slurm job instead of running on the node where `run.sh` was called, which on a
  cluster is a login node. The job's output is in `WORKDIR/slurmout/<stage>_<job>.out`. `EXECUTOR=local` still
  runs them in place. `run.sh all` already ran them inside the finalize job.
- Add `scripts/combine_runs.py` to join validated count tables from separate
  pools, with tests. It refuses unfinished or changed runs, different gene lists
  or count columns, and samples found in more than one run.
- Add a regression test for sample names that begin other names (`MOUSE_1`,
  `MOUSE_10`, `MOUSE_100`). The original concatenation step searched for
  `<sample>*` and kept the first file found, so such a name could receive another
  sample's lane files. Since 2026-09-10 each file is assigned by its full name.
- Point setup documents at the lab's GitHub organization.

## 2026-09-10

- Add sequential Linux execution without Slurm, using the same numbered stages.
- Allow site modules or tools already on PATH; include Ceres and Linux profiles.
- Keep the configured Python selected after loading tool modules.
- Make STAR/index CPU and memory requests configurable for Slurm.
- Record source snapshots, configuration hashes and executable identities.
- Stop later stages if source code or the selected run configuration changes.
- Use the exact input manifest for concatenation, including custom sample names.
- Reject repeated lane directories and overlapping raw/output directories.
- Keep cached reference records when a verification pass is interrupted.
- Validate matching local and Slurm execution on a real two-sample/two-lane subset.
- Clarify setup requirements for the original study-specific R notebooks.

- Add an optional, checksum-verified STAR release setup with a portable local module.
- Record the STAR source commit, binary checksum and license; reuse matching files
  and reject conflicting content. Preserve explicit site-module overrides.
- Include the source MD5 job ID in the submission summary.
- Require matching lane sample sets, complete mates, and unique derived names.
- Record exact input paths and allow an expected sample list.
- Read one FASTQ record for read-length detection.
- Check source MD5 before concatenation in the complete job chain.
- Validate every matrix cell using the selected STAR count column.
- Record per-sample strandedness and reject mixed automatic selections.
- Add HTStream and mapping summaries, BAM validation and indexing to finalization.
- Record validated BAM/index file identities so later changes can be detected.
- Verify reference downloads against published checksums and record index hashes.
- Pin Python dependencies and record installed versions.
- Add separate run settings and recorded Slurm job IDs.
- Preserve count-column, throttle, module and reference-directory overrides.
- Disable scheduler mail by default.
- Update README and setup instructions with inputs, dependencies and acceptance checks.


## 2026-08-26

### Added

- `config.sh`. Every path, name and SLURM setting is here. The same values were
  written into six scripts, so pointing at a
  new study meant editing six files.
- `run.sh`. One entry point for every stage, plus `all` which submits the whole
  pipeline chained with `--dependency=afterok`.
- `run.sh scan`. Reads the lane directories and reports file counts, lane tags,
  read length, and the sample names each naming style would produce. Warns when
  a style would collapse two samples into one name.
- `run.sh check`. Validates the config before anything is submitted: SLURM
  account is real, partition exists, lane directories exist, modules are
  present, workdir is writable. Then works out read length and the STAR index
  path and writes them to `resolved.sh`.
- `References/mouse_rrna.fasta`. 37 KB, required by the HTStream step, not
  previously
  included or described. From the UC Davis Bioinformatics Core
  RNA-Seq workshop the scripts follow.
- `0_Setup/1_fetch_references.slurm`. Downloads the GENCODE genome and
  annotation, verifies the gzip, and checks that the FASTA and GTF use the same
  contig naming before an index is built from them.
- `0_Setup/2_python_env.slurm`. Builds a venv for the two `.py` scripts. No
  module on Ceres carries scipy or pandas.
- `0_Setup/3_verify_source_md5.slurm`. Checks the raw lanes against an
  `@md5Sum.md5` manifest if the sequencing core shipped one.
- `1_Pre_Preprocessing/1_make_samples_txt.sh`. Builds `samples.txt` from the
  complete pairs shared by the lane directories and applies an exclude list.
  The September update changes incomplete lane coverage from an omission to an error.
- `2_HTStream_Preprocessing/3_hts_summary.py`. Reads the per sample HTStream
  JSONs into one table. Nothing read them before.
- `3_STAR_Alignment/6_compare_to_published.py`. Compares a count matrix against
  another one, per sample and per gene.
- `3_STAR_Alignment/7_finalize.slurm`. Runs stats and counts as one job so the
  chain can end on a dependency.
- `EXCLUDE_SAMPLES`. `Undetermined` is excluded by default. The glob picked it
  up as a sample; it is the demultiplexer's unassigned-read bin and is
  usually one of the largest files in the lane.

### Changed

- Concatenation is a SLURM array, one task per sample, instead of one serial
  bash
  loop on the login node. Needs `samples.txt` built first, which is why
  that step moved ahead of it. 173 samples and 970 GB went from serial to 20
  minutes at 35 tasks at once.
- Sample naming is a setting. `cut -d'_' -f1` returns everything before the
  first underscore. On `STUDY_120_S23_L006_R1_001.fastq.gz`
  that is `STUDY` for every sample in the run. `SAMPLE_NAME_STYLE` offers
  `illumina`,
  `first_field` (the previous behavior) and `custom`.
- Lane count and lane tags come from `LANE_DIRS`. Previously fixed at two lanes,
  with `_L008_`
  and `_L004_` written into the find lines.
- `--sjdbOverhang` defaults to read length minus 1, read from the data. It was
  fixed at 100. The index directory is named for the value, so indexes for
  different read lengths can sit side by side.
- The count column is measured instead of assumed. Column 4 of
  `ReadsPerGene.out.tab` was taken
  unconditionally, which fits reverse-stranded libraries only. `stats` sums
  all three columns, picks one, and
  records it.
- SLURM headers carry no account, partition, array size or mail address.
  `run.sh` supplies them on the sbatch command line, where they override the
  file. `#SBATCH` lines are read before the script runs and cannot read a config
  variable, which is why they are passed rather than written.
- `--mail-type` defaults to `END,FAIL`. `ALL` on a 246 task array sends 246
  messages.
- Array size comes from `wc -l samples.txt`. It was hardcoded at
  `--array=1-246`.
- `2_HTStream_Kmer_QC.py` and `5_confirm_reads2counts.py` take arguments and
  read the config from the environment. Both had their paths and sample lists
  written into `main()`.
- `3_star_stats.sh` uses `mktemp -d`. It wrote `names.txt`, `stats.txt`,
  `totals.txt`,
  `temp1` and `temp2` into the current directory and then removed them.
- `4_reads2counts.sh` takes its sample order from `samples.txt` rather than
  `ls`, checks that every sample has the same gene list before pasting, and
  pastes in chunks so it does not run into the open file limit.
- The output matrix filename comes from the config instead of being fixed at
  `GRCm39_AFRI_counts.txt`.
- `2_organize_samples.sh` is gone. Each concatenation task writes straight into
  `concatenated_reads/<sample>/`.
- README restored to the original; separate `SETUP_INSTRUCTIONS.txt` added.

### Fixed

- `1_star_index.slurm` and `2_star.slurm` pointed at different index paths.
  Building the index did not give the alignment step an index it could find.
  Both now read `STAR_INDEX`.
- An empty `$sample` now fails the array task. `sed "${SLURM_ARRAY_TASK_ID}q;d"
  samples.txt`
  returns an empty string when the file is missing, every task then writes to
  `01-HTS_Preproc/.json`, and all of them report success.
- STAR's exit status is checked. The call and its output redirect were inside an
  `eval`'d string, so a crash left an empty BAM and the job exited 0.
- Every stage verifies its own output: gzip integrity and R1/R2 record counts
  after concatenation, expected files after HTStream, a complete index after
  `genomeGenerate`, and column sums against the STAR tables after the matrix is
  built.
- Concatenation checks each selected sample's mates. The sample-list builder
  could still omit samples missing from a lane; this was corrected in September.
- `mkdir` calls that failed on a second run are `mkdir -p`.
- `slurmout/` is created before submitting. Without it SLURM discards the job
  output.

### Notes

- The HTStream chain and the STAR command are unchanged.
- `4_DE_Analysis/*.Rmd` are unchanged. They read an Ensembl annotation export,
  a count table and a study spreadsheet by bare filename; those inputs are not in
  this repository and come from the study's own records.
- The `htstream/1.3.3` module on Ceres contains binaries that report 1.3.2.
- The HTStream link in the previous README returned 404. It is
  https://s4hts.github.io/HTStream/ now.
