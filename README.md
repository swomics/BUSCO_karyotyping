# Quick BUSCO karyotyping :butterfly: :twisted_rightwards_arrows: :butterfly: #
Mapping BUSCO genes from a reference genome to a new assembly to enable more consistent chromosome numbering and crudely infer large scale fissions/fusions

## Dependencies
* Python 3
* (Python 3 module) Biopython
* (Python 3 module) seaborn
* (Python 3 module) reportlab

## Usage

run BUSCO (tested with version 3) for a given reference genome and a query genome. Ensure that chromosomes in the reference are labelled suitably. I used ChrXX for autosomes (e.g. Chr10) and Chr1z for the Z.

For Lepidoptera, I typically run like this:

``` screen -L python run_BUSCO.py -i ./Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -m geno -c 60 -o Hmel2.5 --limit 4 -z -l ~/bin/busco_v3/DB/endopterygota_odb9 -sp heliconius_melpomene1 --long```

### Required options
``` -q (query fasta)```
``` -t (BUSCO "full_table_" output file from query fasta)```
``` -r (BUSCO "full_table_" output file reference fasta)```

### Run with defaults

``` python3.8 BUSCO_2_Chrom.py -q GCA_902806685.1_iAphHyp1.1_genomic.fna -t Query_table.txt -r Hmel2_full_table.txt```


## Outputs :microscope:
Suggested additional chromosome naming based on identity and amount of reference chromosome. For example, LR761654.1_Chr10(36)Chr12(59), denotes the original scaffold name and that it contains 36 BUSCO genes derived from reference chromosome 10 and 59 BUSCO genes derived from chromosome 12.

Additionally, a pdf file is generated displaying a one sided alignment. Specifically, the query genome karyotype along with the specific identity and reference chromosome information of individual BUSCO genes.

![Example output pdf](./output.pdf?raw=true "Example output plot")
