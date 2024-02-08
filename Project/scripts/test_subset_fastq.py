"""Generate the fastq files and map for the test dataset

It subsets 16 randomly chosen Geuvadis samples to reads in the first 2 Mb of chr1.
Those complete samples must already be aligned to be able to subset by region.
"""

import gzip
from pathlib import Path
from subprocess import run

def subset_fastq(fastq: Path, read_ids: set, output: Path):
    """Subset a FASTQ file to a list of read IDs"""
    with gzip.open(fastq, 'rt') as f:
        with gzip.open(output, 'wt') as g:
            for line in f:
                if line.split()[0][1:] in read_ids:
                    g.write(line)
                    for _ in range(3):
                        g.write(next(f))
                else:
                    for _ in range(3):
                        next(f)

sample_file = Path("input/samples.txt")
fqdir = Path("../../Geuvadis/input/fastq")
fqmap = Path("../../Geuvadis/input/fastq_map.txt")
bamdir = Path("../../Geuvadis/intermediate/bam")
outmap = Path("input/fastq_map.txt")
outfqdir = Path("input/fastq")
outfqdir.mkdir(parents=True, exist_ok=True)

with open(sample_file, 'r') as f:
    samples = f.read().splitlines()

with open(fqmap, 'r') as f:
    lines = f.read().splitlines()
with open(outmap, 'w') as out:
    for line in lines:
        fastq1, fastq2, sample = line.split('\t')
        if sample not in samples:
            continue
        print(sample, flush=True)
        bam = bamdir / f'{sample}.bam'
        reads = f'input/fastq/{sample}.chr1_0-2Mb.txt'
        # Get read IDs in chr1:1-2000000 from the bam file
        cmd = f'samtools view {bam} 1:1-2000000 | cut -f1 > {reads}'
        run(cmd, shell=True, check=True)
        with open(reads) as f:
            read_ids = set(f.read().splitlines())
        run(f'rm {reads}', shell=True, check=True)
        # Subset FASTQ file pair to those reads
        out1 = f'chr1_0-2Mb.{fastq1}'
        out2 = f'chr1_0-2Mb.{fastq2}'
        subset_fastq(fqdir / fastq1, read_ids, outfqdir / out1)
        subset_fastq(fqdir / fastq2, read_ids, outfqdir / out2)
        out.write(f'{out1}\t{out2}\t{sample}\n')
