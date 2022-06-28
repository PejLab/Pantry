"""Generate the fastq files and map for the test dataset

It subsets 10 randomly chosen Geuvadis samples to reads in the first Mb of chr1.
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

fqdir = Path("../data/geuvadis/fastq")
fqmap = Path("../data/geuvadis/fastq_map.txt")
bamdir = Path("../out/bam")
outmap = Path("test/fastq_map.txt")
outfqdir = Path("test/fastq")
outfqdir.mkdir(parents=True, exist_ok=True)
# shuf -n 10 data/fastq_map.txt | cut -f3 | sort
samples = [
    'HG00134', 'HG00178', 'HG00267', 'HG00382', 'NA12044',
    'NA18912', 'NA19114', 'NA19223', 'NA20508', 'NA20759',
]

with open(fqmap, 'r') as f:
    lines = f.read().splitlines()
with open(outmap, 'w') as out:
    for line in lines:
        fastq1, fastq2, sample = line.split('\t')
        if sample not in samples:
            continue
        print(sample)
        bam = bamdir / f'{sample}.Aligned.sortedByCoord.out.bam'
        reads = f'test/{sample}.chr1_0-1Mb.txt'
        # Get read IDs in chr1:1-1000000 from the bam file
        cmd = f'samtools view {bam} 1:1-1000000 | cut -f1 > {reads}'
        run(cmd, shell=True, check=True)
        with open(reads) as f:
            read_ids = set(f.read().splitlines())
        # Subset FASTQ file pair to those reads
        out1 = f'chr1_0-1Mb.{fastq1}'
        out2 = f'chr1_0-1Mb.{fastq2}'
        subset_fastq(fqdir / fastq1, read_ids, outfqdir / out1)
        subset_fastq(fqdir / fastq2, read_ids, outfqdir / out2)
        out.write(f'{out1}\t{out2}\t{sample}\n')
