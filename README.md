# ONT_extract
A tool to extract barcode and UMI from Nanopore(Oxford Nanopore Technologies) reads.

# Overview
1. Locate the polyT(or polyA if the read is on the reverse strand) region in the read. Reverse complement the read if it is on the reverse strand. 
2. Discard all the read that is considered as a fused read(contain multiple polyT).
3. Extract cell barcodes and UMIs by using (adapter,linker,polyT) as anchor.
4. Barcode correction. 

# Install
```
pip install celescope editdistance
```

It also needs libssw.so from [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
```
gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h
```
Then put libssw.so under the same directory of `ONT_extract.py`

# Usage
```
usage: ONT_extract.py [-h] [-c CELLBC] [--only_cell] fastq

positional arguments:
  fastq                 fastq file

optional arguments:
  -h, --help            show this help message and exit
  -c CELLBC, --cellBC CELLBC
                        The path of barcodes.tsv file from short-reads library. If provided, the cell barcodes from short reads
                        libary will be used to calculate fraction of reads in cells.
  --only_cell           Only output reads in cells. Use together with --cellBC.
```

# Output
`bc_umi.fa`: fasta file of barcodes and UMIs seperated by ':'
```
>216dc340-3abe-442a-8ac6-d2fac2fb1ff9 88
ATCGACACG_TCTGTCTGC_TGACAGCTA:ATGGCAAGTACT
>fe4f515a-90ea-4ff1-aed4-4ed96ba9f06a 91
GAGCAGCTT_ATAGCGTGT_ACTCGGATT:GAAGTTTTGTTT
>8c84c4a6-787b-4fb5-981f-d86a7aa712c2 84
CGCAACTAC_CCTACAAGG_GAAGCCATT:CCCCGCCAAGTT
>099fae10-de27-4e58-ab70-430c3860db69 88
```

`insert.fq`: fastq file of the insert sequences. Sequences before the polyT are trimmed.

