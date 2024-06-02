#!/usr/bin/python
import os
import sys
import operator
import math
from Bio import SeqIO
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessMultilib(R1file):
  print ("Reading %s" % R1file)
  R2file = R1file.replace('_R1_','_R2_')
  Primerlength = 24
  R1frameshift = 0
  R2frameshift = 0
  roilength    = 105-R1frameshift-R2frameshift
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  variants = [] 
  record_count = 0
  for R1record in R1records:
    record_count += 1
    R2record  = next(R2records)
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi = R1seq[Primerlength+R1frameshift:Primerlength+R1frameshift+roilength]
    R2roi = R2seq[Primerlength+R2frameshift:Primerlength+R2frameshift+roilength]
    if 'N' in R1roi or 'N' in R2roi: continue

    ### We only keep WT from the mutant library
    bc1 = R1roi[84:86] 
    bc2 = rc(R2roi)[84:86]
    if bc1 != bc2: continue
    if bc1 != 'CC' and bc1 != 'TC': continue
    
    ### Only those that have the same forward and reverse reads will be kept
    R1pep = translation(R1roi)
    R2pep = translation(rc(R2roi))
    if R1pep == R2pep:
      variants.append(R1pep) 
    #if record_count == 1000: #break #for small scale testing
      #print (Counter(variants))
      #sys.exit()
  return Counter(variants)

def mut2ID(mut,refseq):
  shift = 25
  haplo = []
  assert(len(mut)==len(refseq))
  for n in range(len(mut)):
    pos = n+shift
    if refseq[n]!=mut[n]:
       haplo.append(refseq[n]+str(pos)+mut[n])
  return '-'.join(haplo)

def Output(inputExp_dict, exp1_dict, exp2_dict, outfile, refseq):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  muts = list(set(list(inputExp_dict.keys())+
                  list(exp1_dict.keys())+
                  list(exp2_dict.keys())))
  WT_exp1_enrich_ratio = float(exp1_dict[refseq]+1)/float(inputExp_dict[refseq]+1)
  WT_exp2_enrich_ratio = float(exp2_dict[refseq]+1)/float(inputExp_dict[refseq]+1)
  outfile.write("\t".join(['mut','mutclass','iptExp_count',
                           'exp1_count','exp2_count',
                           'log_exp1_enrich','log_exp2_enrich',
                           'exp1_fit', 'exp2_fit', 'log_exp1_fit', 'log_exp2_fit'])+"\n")
  for mut in muts:
    mut_exp1_enrich_ratio  = float(exp1_dict[mut]+1)/float(inputExp_dict[mut]+1)
    mut_exp2_enrich_ratio  = float(exp2_dict[mut]+1)/float(inputExp_dict[mut]+1)
    log_mut_exp1_enrich_ratio  = math.log10(float(mut_exp1_enrich_ratio))
    log_mut_exp2_enrich_ratio  = math.log10(float(mut_exp2_enrich_ratio))
    mut_exp1_fitness       = float(mut_exp1_enrich_ratio+1)/float(WT_exp1_enrich_ratio+1)
    mut_exp2_fitness       = float(mut_exp2_enrich_ratio+1)/float(WT_exp2_enrich_ratio+1)
    log_mut_exp1_fitness      = math.log10(float(mut_exp1_fitness))
    log_mut_exp2_fitness      = math.log10(float(mut_exp2_fitness))
    mutclass    = 0 if refseq == mut else hamming(refseq,mut)
    ID          = 'WT' if refseq == mut else mut2ID(mut,refseq)
    outfile.write("\t".join(map(str,[ID,mutclass,inputExp_dict[mut],
                                     exp1_dict[mut],exp2_dict[mut],
                                     log_mut_exp1_enrich_ratio,log_mut_exp2_enrich_ratio,
                                     mut_exp1_fitness,mut_exp2_fitness,log_mut_exp1_fitness,log_mut_exp2_fitness]))+"\n")
  outfile.close()

def main():
  refseq  = next(SeqIO.parse('Fasta/COV107_germline_ref.fa', 'fasta')).seq
  outfile = 'results/COV107_mutlib_fit.tsv'
  inputExp_dict  = ProcessMultilib('fastq/Sample1_ATCACGAT_L001_R1_001.fastq')
  exp1_dict   = ProcessMultilib('fastq/Sample2_CGATGTAT_L001_R1_001.fastq')
  exp2_dict   = ProcessMultilib('fastq/Sample3_TTAGGCAT_L001_R1_001.fastq')
  Output(inputExp_dict, exp1_dict, exp2_dict, outfile, refseq)

if __name__ == "__main__":
  main()
