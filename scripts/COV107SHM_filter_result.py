#!/usr/bin/python
import os
import sys
import string
import operator
import numpy as np
from Bio import SeqIO
from math import log10
from collections import Counter

def convert_ID_seq(mutlist, muts):
  seq = ''
  muts = sorted(list(set([m[0:-1] for m in muts])), key=lambda x:int(x[1::]))
  for mut_ref in muts:
    pos = mut_ref[1::]
    presence = 0
    for mut in mutlist.rsplit('-'):
      if pos in mut:
        seq+=mut[-1]
        presence = 1
        break
    if presence == 0:
      seq+=mut_ref[0]
  return (seq)

def formatting(infile, outfile, inputcutoff, muts):
  #print ("writing: %s" % outfile)
  infile  = open(infile,'r')
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','mutclass','InputExp_Count',
                           'exp1_fit','exp2_fit','exp_avg_fit'])+"\n")
  for line in infile.readlines(): 
    if 'mut' in line: continue
    else:
      line = line.rstrip().rsplit("\t")
      mut  = line[0]
      mutclass = int(line[1])
      iptExp_count = int(line[2])
      exp1_fit   = 10**(float(line[9]))
      exp2_fit   = 10**(float(line[10]))
      if iptExp_count >= inputcutoff: 
        if mut=='WT' or mut.count('-')+1==sum([mut.count(m) for m in muts]):
          exp_avg_fit  = np.mean([exp1_fit,exp2_fit])
          outfile.write("\t".join(map(str,[mut, mutclass, iptExp_count,
                                           exp1_fit, exp2_fit, exp_avg_fit]))+"\n")
          ID = convert_ID_seq(mut, muts)
          print (ID)
  infile.close()
  outfile.close()

def main():
  inputcutoff = 0
  muts    = ['G26E','F27L','F27V','F27I','T28I','S31R','S35T','V50L','S53P','S56T','T57A','Y58F']
  infile  = 'results/COV107_mutlib_fit.tsv'
  outfile = 'results/COV107_mutlib_fit_exp.tsv'
  formatting(infile, outfile, inputcutoff, muts)

if __name__ == "__main__":
  main()
