import multiprocessing as mp
import sys
import glob
from functools import partial
import subprocess

def generate_chromosome_list(regions_path):
    return [i.rsplit("/", 1)[1].replace(".regions.txt", "") for i in glob.glob(f"{regions_path}/*.regions.txt")]

def load_lengths(filename="/mnt/staging/hg19/hg19.fa.len"):
    lens = {}
    with open(filename, 'r') as f:
        for row in f:
            chrom, length = row.strip().split("\t")
            lens[chrom] = int(length)
    return lens

def process_chromosome(chrom, max_length_dict, bam_file, sample_name, step_size=250000):
    HEADER=True
    max_length=max_length_dict.get(chrom, None)
    if max_length:
        for i in range(1, max_length, step_size):
            if i+(step_size-1) < max_length:
                if HEADER:
                    CMD = f"/data/cgreco/bin/samtools coverage -r {chrom}:{i}-{i+(step_size-1)} {bam_file} >> {sample_name}.{chrom}.cov"
                    HEADER=False
                else:
                    CMD = f"/data/cgreco/bin/samtools coverage -H -r {chrom}:{i}-{i+(step_size-1)} {bam_file} >> {sample_name}.{chrom}.cov"
            else:
                if HEADER:
                    CMD = f"/data/cgreco/bin/samtools coverage -r {chrom}:{i} {bam_file} >> {sample_name}.{chrom}.cov"
                    HEADER=False
                else:
                    CMD = f"/data/cgreco/bin/samtools coverage -H -r {chrom}:{i} {bam_file} >> {sample_name}.{chrom}.cov"
            subprocess.call(CMD, shell=True)
    else:
        # Captures non-canonical (_random, hap, etc.) fasta's
        print({chrom})
        CMD = f"/data/cgreco/bin/samtools coverage -r {chrom} {bam_file} >> {sample_name}.{chrom}.cov"
        subprocess.call(CMD,shell=True)
            

def main():
    regions_dir = sys.argv[1]
    bam_file = sys.argv[2]
    threads = int(sys.argv[3])
    sample_name = sys.argv[4]
    chrom_list = generate_chromosome_list(regions_dir)
    max_length_dict = load_lengths()
    with mp.Pool(threads) as p:
        p.map(partial(process_chromosome, max_length_dict=max_length_dict, bam_file=bam_file, sample_name=sample_name), chrom_list)

if __name__=="__main__":
    main()


