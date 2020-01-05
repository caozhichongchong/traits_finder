# traits_finder
## Introduction
* traits_finder searches and summarizes traits in genomes and metagenomes
* input: reference database and folder of genomes/metagenomes
* requirement: python >= 3.0, blast or hmm
* requirement: for hmm, you need to prepare the hmm database
* optional: diamond, bwa, hs-blastn, usearch, mafft, fasttree

## Install
`pip install traits_finder`\
in preparation: `anaconda download caozhichongchong/traits_finder`
### latest version (unstable though)
`git clone https://github.com/caozhichongchong/traits_finder.git `\
`cd traits_finder`\
`python setup.py build`\
`python setup.py install`

## Availability

https://pypi.org/project/traits_finder

## What do you need to prepare
1. your reference database (-db your.db), protein sequences (-dbf 1) or dna sequences (-dbf 2)
2. a mapping file of functions to each reference sequence (sequence	function)
3. all genomes/metagenomes in a folder (-i your.input.folder)
4. suffix or file extension of your genomes/metagenomes, such as .fasta or .fastq (-fa your.input.genome/metagenome.format)
5. programs to run: blast for similarity search (-s 1 )or hmm for domain search (-s 2)
6. optional programs to speedup! diamond, hs-blastn (or usearch for 16S extracting), usearch (necessary for HGT and HGT_sum)
7. optional programs to look at sequence variants! bwa, mafft, fasttree

## How to use it
### Search traits: boring and slow...
#### Database installed in traits_finder: antibiotic resistant genes (-db ARG), butyrate producing genes (-db but)
1. search protein reference sequences in genomes (traits_finder genome) or mobile genetic elements (traits_finder mge) by similarity search\
`traits_finder genome -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u diamond --bp blastp -dbf 1 -s 1`\

2. search protein reference sequences in metagenomes by similarity search\
`traits_finder meta -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u diamond --bp blastp -dbf 1 -s 1`\

3. search dna reference sequences in genomes by similarity search\
`traits_finder genome -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn -dbf 2 -s 1`\

4. search dba reference sequences in metagenomes by similarity search\
`traits_finder meta -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn -dbf 2 -s 1`\

5. search protein reference sequences in genomes by hmm\
`traits_finder genome -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --hmm hmmscan -dbf 1 -s 2`\

6. search dna reference sequences in genomes by alignment\
`traits_finder genome -db your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn --bwa bwa -dbf 2 -s 1`\

7. search dna reference sequences in metagenomes by alignment\
`traits_finder meta -db your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --u usearch.or.hs-blastn --bp blastn --bwa bwa -dbf 2 -s 1`\

### Summarize results: cool and fast!
1. summarize traits in genome\
`traits_finder sum_genome -db your.db -m function.mapping.your.db -i your.input.folder -fa your.input.genome.format --orf your.input.orf.format --r your.output.folder --r16 your.output.folder.for.16s`\

2. summarize traits in metagenomes\
`traits_finder sum_meta -db your.db -m function.mapping.your.db -i your.input.folder -fa your.input.metagenomes.format --r your.output.folder --r16 your.output.folder.for.16s --meta metadata.txt`\

### HGT finder and summarizing: cool and fast! (still-testing)

## Results

## Copyright
Copyright: An Ni Zhang, Prof. Eric Alm, Alm Lab in MIT\
Citation: Not yet, coming soon!\
Contact: anniz44@mit.edu
