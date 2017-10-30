#!/usr/bin/env python3

import re
import os
import gffutils
from sys import argv
import subprocess
import gffutils.gffwriter as gffwriter
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from multiprocessing import Pool
from Bio import Seq
import regex as re


def overlap(gene, hint):
    gff_gene = open(gene, "r")
    gff_hint = open(hint, "r")
    hints = {}
    genes = {}
    print ("sample overlap genes hints sensitivitys specificity score")

    
    for line in gff_hint:
        fields = line.split("\t")
        if int(fields[5]) > 2:
            uniq = fields[0] + fields[3] +fields[4] + fields[6]
            hints[uniq] = line
    for line in gff_gene:
        fields = line.split("\t")
        uniq = fields[0] + fields[3] +fields[4] + fields[6]
        genes[uniq] = line
    overlap = {}
    for key in hints:
        if key in genes:
            overlap[key] = genes[key]
    print (gene , len(overlap), len(genes), len(hints), len(overlap)/len(genes)*100 , len(overlap)/len(hints)*100, (len(overlap)/len(genes)*100 + len(overlap)/len(hints)*100) /2)

if __name__ == '__main__':
    
    
    gene = argv[1]
    hint = argv[2] 
    overlap(gene, hint)
    