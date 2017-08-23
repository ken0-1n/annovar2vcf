#! /usr/bin/env/python

import sys
import vcf
from pyfasta import Fasta

mutations_anno = sys.argv[1]
f = Fasta('../reference/GRCh37.fa')
vcf_reader = vcf.Reader(filename='./template.vcf')
vcf_writer = vcf.Writer(open('./out.vcf', 'w'),vcf_reader)

header_flag = True
with open(mutations_anno, 'r') as hin:
    for line in hin:
        if line.startswith("#"): continue
        if header_flag:
            header_flag = False
            continue
        line = line.rstrip('\n')
        F = line.split('\t')
        chrom = F[0]
        start = F[1]
        end = F[2]
        ref = F[3]
        alt = F[4]

        fchrom = chrom

        if ref == "-":
            # print chrom+"\t"+start+"\t"+end+"\t"+ref+"\t"+alt
            fstart = int(start)-1
            fend = int(end)
            ref_seq = f[fchrom][fstart:fend]
            # print chrom+"\t"+start+"\t"+end+"\t"+ref_seq+"\t"+ref_seq+alt

            sub_record = vcf.model._Substitution(ref_seq+alt)
            record = vcf.model._Record(chrom,int(start),None,ref_seq,[sub_record],None,[],None,None,None)

            vcf_writer.write_record(record)

        elif alt == "-":
            # print chrom+"\t"+start+"\t"+end+"\t"+ref+"\t"+alt
            fstart = int(start)-2
            fend = int(start)-1
            ref_seq = f[fchrom][fstart:fend]
            # print chrom+"\t"+str(int(start)-1)+"\t"+end+"\t"+ref_seq+ref+"\t"+ref_seq

            sub_record = vcf.model._Substitution(ref_seq)
            record = vcf.model._Record(chrom,int(start)-1,None,ref_seq+ref,[sub_record],None,[],None,None,None)

            vcf_writer.write_record(record)

vcf_writer.close()        


