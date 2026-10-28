import gzip,importlib.util,json,os,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(path):
    spec=importlib.util.spec_from_file_location(path.stem,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
sample=load(ROOT/'scripts/sample_manifest.py');counts=load(ROOT/'3_STAR_Alignment/5_confirm_reads2counts.py')
class Validation(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
    def lanes(self):
        lanes=[]
        for lane in ['l1','l2']:
            d=self.root/lane;d.mkdir();lanes.append(d)
            for s in ['STUDY_2','STUDY_10','Undetermined']:
                for mate in [1,2]:
                    with gzip.open(d/f'{s}_S1_L001_R{mate}_001.fastq.gz','wt') as f:f.write('@a\nACGT\n+\nIIII\n')
        return lanes
    def collect(self,lanes):
        return sample.collect(lanes,'s/_S[0-9]+_L[0-9]+$//',['_R1_001.fastq.gz','_R2_001.fastq.gz'],{'Undetermined'})
    def test_two_lanes(self):
        names,rows=self.collect(self.lanes());self.assertEqual(names,['STUDY_2','STUDY_10']);self.assertEqual(len(rows),4)
    def test_missing_mate(self):
        lanes=self.lanes();next(lanes[1].glob('STUDY_2*R2*')).unlink()
        with self.assertRaises(ValueError):self.collect(lanes)
    def test_missing_whole_sample(self):
        lanes=self.lanes()
        for p in lanes[1].glob('STUDY_2*'):p.unlink()
        with self.assertRaises(ValueError):self.collect(lanes)
    def test_collision(self):
        lanes=self.lanes();p=lanes[0]/'STUDY_2_S2_L001_R1_001.fastq.gz';p.write_bytes(b'x')
        with self.assertRaises(ValueError):self.collect(lanes)
    def test_prefix_names_keep_their_own_files(self):
        # A name that starts another name (S_1, S_10, S_100) must still get only its own files.
        lanes=[]
        for lane,number in [('l1',5),('l2',6)]:
            d=self.root/lane;d.mkdir();lanes.append(d)
            for i,s in enumerate(['MOUSE_1','MOUSE_10','MOUSE_100','MOUSE_101'],1):
                for mate in [1,2]:
                    with gzip.open(d/f'{s}_S{i}_L00{number}_R{mate}_001.fastq.gz','wt') as f:f.write(f'@{s}\nACGT\n+\nIIII\n')
        names,rows=self.collect(lanes)
        self.assertEqual(names,['MOUSE_1','MOUSE_10','MOUSE_100','MOUSE_101'])
        self.assertEqual(len(rows),8)
        for name,lane,r1,r2,*_ in rows:
            self.assertTrue(Path(r1).name.startswith(name+'_S'));self.assertTrue(Path(r2).name.startswith(name+'_S'))
    def test_read_length_first_record(self):
        p=self.root/'reads.gz'
        with gzip.open(p,'wt') as f:f.write('@a\nACGT\n+\nIIII\n'+'not another FASTQ record\n'*10000)
        self.assertEqual(sample.read_length(p),4)
    def fixtures(self,column):
        star=self.root/'star';star.mkdir();s=star/'s';s.mkdir();samples=self.root/'samples';samples.write_text('s\n')
        (s/'s_ReadsPerGene.out.tab').write_text(''.join(x+'\t0\t0\t0\n' for x in ['N_unmapped','N_multimapping','N_noFeature','N_ambiguous'])+'gene1\t3\t4\t5\ngene2\t6\t7\t8\n')
        matrix=self.root/'counts';matrix.write_text(f'gene_id\ts\ngene1\t{column+1}\ngene2\t{column+4}\n')
        return star,matrix,samples
    def test_columns(self):
        for column in [2,3,4]:
            with self.subTest(column=column):
                if (self.root/'star').exists():
                    import shutil;shutil.rmtree(self.root/'star')
                args=self.fixtures(column);self.assertEqual(counts.compare_counts(*args,column)['cells_compared'],2)
    def test_equal_sum_swapped_counts(self):
        args=self.fixtures(4);args[1].write_text('gene_id\ts\ngene1\t8\ngene2\t5\n')
        with self.assertRaises(ValueError):counts.compare_counts(*args,4)
    def test_missing_sample(self):
        args=self.fixtures(4);args[2].write_text('s\nmissing\n')
        with self.assertRaises(ValueError):counts.compare_counts(*args,4)
    def test_duplicate_gene(self):
        args=self.fixtures(4);args[1].write_text('gene_id\ts\ngene1\t5\ngene1\t8\n')
        with self.assertRaises(ValueError):counts.compare_counts(*args,4)
    def test_run_settings_override_defaults(self):
        settings=self.root/'settings.sh'
        settings.write_text(f'WORKDIR="{self.root}/run"\nCOUNT_COLUMN=3\nTHROTTLE_STAR=2\nMODULE_STAR="star/test"\n')
        env=os.environ.copy();env['PIPELINE_CONFIG']=str(settings)
        command='source "$1"; printf "%s\\n" "$REF_DIR" "$COUNT_COLUMN" "$THROTTLE_STAR" "$MODULE_STAR"'
        result=subprocess.check_output(['bash','-c',command,'bash',str(ROOT/'config.sh')],env=env,text=True)
        self.assertEqual(result.splitlines(),[str(self.root/'run/References'),'3','2','star/test'])
        settings.write_text(settings.read_text()+f'REF_DIR="{self.root}/shared_refs"\n')
        result=subprocess.check_output(['bash','-c',command,'bash',str(ROOT/'config.sh')],env=env,text=True)
        self.assertEqual(result.splitlines()[0],str(self.root/'shared_refs'))
    def test_final_qc_records_validated_file_identity(self):
        (self.root/'samples.txt').write_text('sample\n')
        align=self.root/'02-STAR_alignment/sample';align.mkdir(parents=True)
        bam=align/'sample_Aligned.sortedByCoord.out.bam';bam.write_bytes(b'fixture')
        fields=['Number of input reads','Uniquely mapped reads number','Uniquely mapped reads %','Number of reads mapped to multiple loci','% of reads mapped to multiple loci','% of reads unmapped: too short']
        (align/'sample_Log.final.out').write_text('\n'.join(x+' | 10' for x in fields))
        bin_dir=self.root/'bin';bin_dir.mkdir();stub=bin_dir/'samtools'
        stub.write_text('#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nif sys.argv[1]=="index":Path(sys.argv[-1]+".bai").write_bytes(b"index fixture")\n');stub.chmod(0o700)
        env=os.environ.copy();env['PATH']=str(bin_dir)+os.pathsep+env['PATH']
        subprocess.run(['python3',str(ROOT/'scripts/final_qc.py'),str(self.root)],env=env,check=True,capture_output=True)
        report=json.loads((self.root/'bam_validation.json').read_text());self.assertEqual(len(report['files']),2)
        for name,value in report['files'].items():
            stat=(self.root/name).stat();self.assertEqual(value,[stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns,stat.st_ino])
if __name__=='__main__':unittest.main()
