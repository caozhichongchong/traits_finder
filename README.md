# traits_finder
## Introduction
* traits_finder searches and summarizes traits in genomes and metagenomes
* input: reference database and folder of genomes/metagenomes
* requirement: blast or hmm
* requirement: for hmm, you need to prepare the hmm database
* optional: diamond, bwa, hs-blastn, usearch

## Install
`pip install traits_finder`\
in preparation: `anaconda download caozhichongchong/traits_finder`

## Availability

https://pypi.org/project/traits_finder

## How to use it
1. search protein reference sequences in genomes by similarity search\
`traits_finder -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u diamond --bp blastp -dbf 1 -inf 1 -s 1`\

2. search protein reference sequences in metagenomes by similarity search\
`traits_finder -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u diamond --bp blastp -dbf 1 -inf 2 -s 1`\

3. search dna reference sequences in genomes by similarity search\
`traits_finder -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn -dbf 2 -inf 1 -s 1`\

4. search dba reference sequences in metagenomes by similarity search\
`traits_finder -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn -dbf 2 -inf 2 -s 1`\

5. search protein reference sequences in genomes by hmm\
`traits_finder -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --hmm hmmscan -dbf 1 -inf 1 -s 2`\

6. search dna reference sequences in genomes by alignment\
`traits_finder -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn --bwa bwa -dbf 2 -inf 1 -s 1`\

7. search dna reference sequences in metagenomes by alignment\
`traits_finder -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn --bwa bwa -dbf 2 -inf 2 -s 1`\

## Results

## Copyright
Copyright: An Ni Zhang, Prof. Eric Alm, Alm Lab in MIT\
Citation: Not yet, coming soon!\
Contact: anniz44@mit.edu
