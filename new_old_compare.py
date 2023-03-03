#!/usr/bin/env python 

import argparse
import pysam
import vcf
import pyfaidx
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--vcf")
parser.add_argument("--ground_truth")
parser.add_argument("--fasta")
parser.add_argument("--output")
parser.add_argument("--bam")
args = parser.parse_args()

fasta = pyfaidx.Fasta(args.fasta)
vcfin = vcf.Reader(filename =args.vcf)
ground_truth = vcf.Reader(filename=args.ground_truth)
bam = pysam.AlignmentFile(args.bam)

def pair_iter(i1, i2, key):
    ''' Iterate through sorted i1 and i2 simulateneously. Emit (v1, v2) for items
        that compare equal. Emit (v1, None) for value v1 from i1 that doesn't have
        an equal item in i2. Emit (None, v2) for value v2 from i2 that doesn't have
        and equal items in i1.  Comparisons are made with key(v).  Assumes
        that i1 and i2 are sorted with respect to key(). Also assumes that
        a sigle iterator does not contain items equal to each other.'''

    v1 = None
    v2 = None

    while True:

        if v1 is None:
            try:
                v1 = next(i1)#i1.next()
            except StopIteration:
                # return the rest of i2
                if v2 is not None:
                    yield (None, v2)
                for x2 in i2:
                    yield (None, x2)
                break

        if v2 is None:
            try:
                v2 = next(i2) #i2.next()
            except StopIteration:
                # return the rest of i2
                if v1 is not None:
                    yield (v1, None)
                for x1 in i1:
                    yield (x1, None)
                break

        k1 = key(v1)
        k2 = key(v2)
        if k1 == k2:
            yield (v1, v2)
            v1 = None
            v2 = None
        elif k1 < k2:
            yield (v1, None)
            v1 = None
        else:
            yield (None, v2)
            v2 = None

output = open(args.output, 'w')
output.write(",".join(["chrom","pos", "in_vcf","in_gt","qual","ref","alt","is_indel","depth","AD", "PS_vcf", "PS_gt", "CC1", "CC2", "MS1", "MS2"])+"\n")
for chrom in fasta.keys():
    try:
        vcfin.fetch(chrom)
    except:
        print("could not fetch chrom in vcf",chrom)
        continue
    try:
        ground_truth.fetch(chrom)
    except:
        print("failed to fetch in gt",chrom)
        continue
    
    
    for (vcf_rec, gt_rec) in pair_iter(vcfin, ground_truth, lambda x: x.POS):
        print(vcf_rec, gt_rec)
        pos = "0"
        in_vcf = "0"
        in_gt = "0"
        qual = "0"
        ref = ""
        alt = ""
        is_indel = "0"
        depth = "-1"
        AD = "-1"
        phase_set = "-1"
        cluster_center1 = "-1"
        cluster_center2 = "-1"
        molecule_support1 = "-1"
        molecule_support2 = "-1"
        if vcf_rec:
            pos = str(vcf_rec.POS)
            in_vcf = "1"
            qual = str(vcf_rec.QUAL)
            ref = str(vcf_rec.REF)
            alt = str(vcf_rec.ALT)
            if vcf_rec.is_indel():
                is_indel = "1"
            depth = str(vcf_rec.samples[0]["DP"])
            AD = str(vcf_rec.samples[0]["AD"][1]) 
            if "PS" in vcf_rec.samples[0]:
                phase_set = str(vcf_rec.samples[0]["PS"])
            if "CC" in vcf_rec.samples[0]:
                toks = vcf_rec.samples[0]["CC"]
                cluster_center1 = str(toks[0])
                cluster_center2 = str(toks[1])
            if "MS" in vcf_rec.samples[0]:
                toks = vcf_rec.samples[0]["MS"]
                molecule_support1 = str(toks[0])
                molecule_support2 = str(toks[1])
        else:
            depth = 0
            for read in bam.fetch(chrom, pos, pos+1):
                depth += 1
            depth = str(depth)
            
        if gt_rec:
            pos = str(gt_rec.POS)
            in_gt = "1"
            ref = str(vcf_rec.REF)
            alt = str(vcf_rec.ALT[0])
            gt_phase_set = "-1"
            if gt_rec.is_indel():
                is_indel = "1"
            if "PS" in gt_rec.samples[0]:
                gt_phase_set = gt_rec.samples[0]["PS"]


        output.write(",".join([chrom,pos,in_vcf,in_gt,qual,ref,alt,is_indel,depth,AD,phase_set,gt_phase_set,
            cluster_center1, cluster_center2, molecule_support1, molecule_support2])+"\n")
