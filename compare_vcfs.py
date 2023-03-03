#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='compare vcfs and output paired listings with depth and other info.')
parser.add_argument('--vcf')
parser.add_argument('--ground_truth')
parser.add_argument('--output')
parser.add_argument('--bam')
parser.add_argument('--regions') 
parser.add_argument('--fasta')

args = parser.parse_args()

import vcf
import pyfaidx
import pysam
print(args.vcf)
print(args.ground_truth)
#bam = pysam.AlignmentFile(args.bam)

fasta = pyfaidx.Fasta(args.fasta)

def pair_iter(i1, i2):
    v1 = None
    v2 = None

    while True:
        if v1 is None:
            try:
                v1 = next(i1)
            except StopIteration:
                if v2 is not None:
                    yield (None, v2)
                for x2 in i2:
                    yield (None, x2)
                break
        if v2 is None:
            try:
                v2 = next(i2)
            except StopIteration:
                if v1 is not None:
                    yield (v1, None)
                for x1 in i1:
                    yield (x1, None)
                break
        k1 = v1.POS
        k2 = v2.POS
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



vcf1 = vcf.Reader(filename=args.vcf)
vcf2 = vcf.Reader(filename=args.ground_truth)
def get_gt(rec, sample):
    gt_bases = sample.gt_bases
    gt_bases = gt_bases.split("/")
    gt_bases = gt_bases.split("|")
    ref = rec.alleles[0]
    separator = "/"
    if sample.phased:
        separator = "|"
    bases = []
    for x in gt_bases:
        if x == ref:
            bases.append("0")
        else:
            bases.append("1")
    gt_vcf = separator.join(bases)
    return(gt_vcf)

import numpy as np
from sklearn.neighbors import NearestNeighbors

closest_indel = {}
'''
for contig in fasta.keys():
    print("building nearest neighbor for "+str(contig))
    try:
        vcf1.fetch(contig)
    except:
        continue
    indel_positions = []
    for rec in vcf1:
        if rec.is_indel:
            indel_positions.append([rec.POS]) 
    neighbors = NearestNeighbors(n_neighbors=1)    
    if len(indel_positions) > 0:
        neighbors.fit(indel_positions)
        closest_indel[contig] = neighbors
'''
with open(args.output,'w') as out:
    out.write("\t".join(["chrom","pos","ref","alt","is_indel","phased","closest_indel","in_vcf","in_gt","depth","is_het_vcf","is_het_gt","gt_vcf","gt_gt","ps_vcf","ps_gt","hap1","hap2","hap1mols","hap2mols"])+"\n")
    for contig in fasta.keys():
        try:
            print("fetching "+str(contig))
            vcf1.fetch(contig)
            vcf2.fetch(contig)
        except:
            print("failed fetching "+str(contig))
            continue
        for (putative, ground_truth) in pair_iter(vcf1, vcf2):
            is_het_vcf = "0"
            is_het_gt = "0"
            in_vcf = "0"
            in_gt = "0"
            pos = 0
            gt_vcf = "./."
            gt_gt = "./."
            ps_vcf = "-1"
            ps_gt = "-1"
            hap1 = "-1"
            hap2 = "-1"
            hap1mols = "-1"
            hap2mols = "-1"
            phased = "0"
            ref = "N"
            alt = "N"
            is_indel = "0"
            indel_dist = "-1"
            if putative:
                in_vcf = "1"
                pos = putative.POS
                if putative.CHROM in closest_indel:
                    indel_dist = closest_indel[putative.CHROM].kneighbors([[pos]],1,return_distance=True)[0][0][0] # I have no idea. I hope it doesnt break
                ref = putative.REF
                alt = putative.ALT[0]
                if putative.is_indel:
                    is_indel = "1"
                
                sample = putative.samples[0]
                if sample.is_het:
                    is_het_vcf = "1"
                gt_vcf = sample['GT']#get_gt(putative, sample)
                if sample.phased:
                    phased = "1"
                try:
                    haps = sample['CC']
                    hap1 = str(haps[0])
                    hap2 = str(haps[1])
                except:
                    pass
                try:
                    mols = sample["MS"]
                    hap1mols = str(mols[0])
                    hap2mols = str(mols[1])
                except:
                    pass
                try:
                    ps_vcf = str(sample["PS"])
                except:
                    pass
            if ground_truth:
                in_gt = "1"
                pos = ground_truth.POS
                ref = ground_truth.REF
                alt = ground_truth.ALT[0]
                if ground_truth.is_indel:
                    is_indel = "1"
                if ground_truth.CHROM in closest_indel:
                    indel_dist = closest_indel[ground_truth.CHROM].kneighbors([[pos]],1,return_distance=True)[0][0][0]
                sample = ground_truth.samples[0]
                if sample.is_het:
                    is_het_gt = "1"
                gt_gt = sample['GT']#get_gt(ground_truth, sample)
                try:
                    ps_gt = str(sample["PS"])
                except:
                    pass
            depth = 0
            #for read in bam.fetch(contig, pos, pos+1):
            #    depth += 1
            
            out.write("\t".join([contig, str(pos), str(ref), str(alt), phased, is_indel, str(indel_dist), in_vcf, 
                in_gt, str(depth), is_het_vcf, is_het_gt,
                gt_vcf, gt_gt, ps_vcf, ps_gt, hap1, hap2, hap1mols, hap2mols])+"\n")
            out.write("\t".join(["chrom","pos","ref","alt","is_indel","phased","closest_indel","in_vcf","in_gt","depth","is_het_vcf","is_het_gt","gt_vcf","gt_gt","ps_vcf","ps_gt","hap1","hap2","hap1mols","hap2mols"])+"\n")
