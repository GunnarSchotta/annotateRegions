#!/usr/bin/env python3
"""Generate a random BED file of test regions on the mm10 assembly."""

import random

random.seed(42)

MM10_CHROM_SIZES = {
    "chr1": 195471971,
    "chr2": 182113224,
    "chr3": 160039680,
    "chr4": 156508116,
    "chr5": 151834684,
    "chr6": 149736546,
    "chr7": 145441459,
    "chr8": 129401213,
    "chr9": 124595110,
    "chr10": 130694993,
    "chr11": 122082543,
    "chr12": 120129022,
    "chr13": 120421639,
    "chr14": 124902244,
    "chr15": 104043685,
    "chr16": 98207768,
    "chr17": 94987271,
    "chr18": 90702639,
    "chr19": 61431566,
    "chrX": 171031299,
}

N_REGIONS = 200
MIN_WIDTH = 200
MAX_WIDTH = 800

chroms = list(MM10_CHROM_SIZES.keys())
regions = []
for i in range(N_REGIONS):
    chrom = random.choice(chroms)
    width = random.randint(MIN_WIDTH, MAX_WIDTH)
    start = random.randint(0, MM10_CHROM_SIZES[chrom] - width - 1)
    end = start + width
    regions.append((chrom, start, end, f"region_{i+1}"))

regions.sort(key=lambda r: (chroms.index(r[0]), r[1]))

with open("regions.bed", "w") as f:
    for chrom, start, end, name in regions:
        f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

print(f"Wrote {len(regions)} regions to regions.bed")
