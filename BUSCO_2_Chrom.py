#!/usr/bin/env python3
#Author: Samuel Whiteford

"""
This script maps BUSCOs coloured by reference (ancestral) chromosome to an input assembly. The output is a karyotype plot and suggested new scaffold/chromosome names.
"""


import argparse
import sys
from Bio import SeqIO
from Bio.Graphics import BasicChromosome
from reportlab.lib.units import cm
from Bio.SeqFeature import SeqFeature, FeatureLocation
import seaborn as sns
from collections import Counter

parser = argparse.ArgumentParser()

parser.add_argument('-r', action='store',
                    dest='reftable',
                    help='Path to reference table')

parser.add_argument('-q', action='store',
                    dest='queryfasta',
                    help='Path to query fasta')

parser.add_argument('-t', action='store',
                    dest='querytable',
                    help='Path to query table')

parser.add_argument('-o', action='store',
                    dest='out', default="output.pdf",
                    help='Minimum length to consider a sequence as a chromosome')

parser.add_argument('-m', action='store',
                    dest='min', default=1000000,
                    help='Minimum length to consider a sequence as a chromosome')

parser.add_argument('-b', action='store',
                    dest='minb', default=10,
                    help='Minimum number of BUSCOs to assign ancestry')

parser.add_argument('--verbose',dest='verb', action='store_true', help='print additional output')

parser.add_argument('--version', action='version',
                    version='%(prog)s 1.0')

args = parser.parse_args()

## function to filter dictionary
def filterTheDict(dictObj, callback):
    newDict = dict()
    # Iterate over all the items in dictionary
    for (key, value) in dictObj.items():
        # Check if item satisfies the given condition then add to new dict
        if callback((key, value)):
            newDict[key] = value
    return newDict


# initialise max length variable and empty dictionaries for storing information on input fasta, reference busco table and query busco table
max_len=0
karyotype_dict = dict()
karyotype_content_dict = dict()
BUSCO_ref_dict = dict()
BUSCO_query_dict = dict()
colour_dict = dict()

# initialise karyotype plot variables
telomere_length = 200000
chr_diagram = BasicChromosome.Organism()
chr_diagram.page_size = (120 * cm, 60 * cm)  # A4 landscape

##initialise features and feature dict
out=[]
feat_dict = dict()

print("Parsing fasta",end="\t")
### parse query fasta file into names and lengths
### set max length
for record in SeqIO.parse(args.queryfasta,"fasta"):
    if len(record.seq) > max_len:
        max_len = len(record.seq)
    if len(record.seq) > args.min:
        karyotype_content_dict[record.id] = []
        karyotype_dict[record.id] = len(record.seq)
        if args.verb==True: print(record.id, record.description, "length:", len(record.seq))
print("Done")
print(len(karyotype_dict),"sequences over the minimum size threshold (-m)",args.min)

## Do some basic checks on query fasta
if bool(karyotype_dict)==False:
    sys.exit("ERROR: No sequences above minimum length threshold")
if args.min < 0 and isinstance(args.min, int)==False:
    sys.exit("ERROR: Minimum size threshold must be a positive integer")

## Parse reference BUSCO table and initialise colours
print("Parsing BUSCO tables",end="\t")
BUSCO_ref = open(args.reftable)
for line in BUSCO_ref:
    if line.startswith('#'):
        continue
    else:
        line2 = line.rstrip()
        line_array = line2.split("\t")
        if line_array[1] == "Complete":
            BUSCO_ref_dict[line_array[0]] = line_array[2]
            colour_dict[line_array[2]] = 0
            if args.verb == True: print(line_array)

## TBD do some basic checks on reference BUSCO e.g. high number of missing fragmented etc

## set colours for reference LGs
palette = sns.color_palette("hls", len(colour_dict))
for LG in colour_dict:
    colour_dict[LG] = palette.pop(0)

## Parse query BUSCO table and write features
BUSCO_query = open(args.querytable)
for line in BUSCO_query:
    if line.startswith('#'):
        continue
    else:
        line2 = line.rstrip()
        line_array = line2.split("\t")
        if line_array[1] == "Complete" and line_array[0] in BUSCO_ref_dict and line_array[2] in karyotype_dict:
            BUSCO_query_dict[line_array[0]] = line_array
            feat_dict[line_array[2]]=[]
            karyotype_content_dict[line_array[2]].append(BUSCO_ref_dict[line_array[0]])
            if args.verb == True: print(line_array)
print("Done")

## create feature set for mapped BUSCOs
print("Building one-sided alignment karyotype plot")
for qBusco in BUSCO_query_dict:
    arr = BUSCO_query_dict[qBusco]
    LG = BUSCO_ref_dict[arr[0]]
    feature = SeqFeature(FeatureLocation(int(arr[3]), int(arr[4]), ref=arr[2]), strand=0,
                         qualifiers={"locus_tag": [LG+"_"+arr[0]], "color": [colour_dict[LG]]})
    feat_dict[arr[2]].append(feature)
    #out.append(feature)

# for record in
#     record = SeqIO.read(filename, "fasta")
#     print(name, len(record))
#
# max_len = 30432563  # Could compute this from the entries dict
# For illustration

#features = dict()
# for index, (name, filename) in enumerate(entries):
for key in karyotype_dict:

    name = key
    length = karyotype_dict[key]

    #     record = SeqIO.read(filename, "genbank")
    #     length = len(record)
    #     features = [f for f in record.features if f.type == "tRNA"]
    #     # Record an Artemis style integer color in the feature's qualifiers,
    #     # 1 = Black, 2 = Red, 3 = Green, 4 = blue, 5 =cyan, 6 = purple
    #     for f in features:
    #         f.qualifiers["color"] = [index + 2]
    #
    cur_chromosome = BasicChromosome.Chromosome(name)
    #     # Set the scale to the MAXIMUM length plus the two telomeres in bp,
    #     # want the same scale used on all five chromosomes so they can be
    #     # compared to each other
    cur_chromosome.scale_num = max_len + 2 * telomere_length
    #
    #     # Add an opening telomere
    start = BasicChromosome.TelomereSegment()
    start.scale = telomere_length
    cur_chromosome.add(start)
    #
    #     # Add a body - again using bp as the scale length here.
    #     body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
    body = BasicChromosome.AnnotatedChromosomeSegment(length,feat_dict[name])
    body.scale = length
    cur_chromosome.add(body)
    #
    #     # Add a closing telomere
    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_length
    cur_chromosome.add(end)
    #
    #     # This chromosome is done
    chr_diagram.add(cur_chromosome)
#
chr_diagram.draw(args.out, args.queryfasta)
#


for key in karyotype_content_dict:
    count = Counter(karyotype_content_dict[key])
    newDict = filterTheDict(count, lambda elem: elem[1] >= args.minb)
    print(key,key,end="_",sep="\t")
    for item, amount in newDict.items(): 
        print("{}({})".format(item, amount),end="")
    print("\t",key,sep="\t",end="_")
    for item, amount in newDict.items(): 
        print("{}".format(item),end="")
    print()
