#!/usr/bin/env python

import argparse
import pandas as pd
import sys
parser = argparse.ArgumentParser()
parser.add_argument("--output")
parser.add_argument("--csv")
args = parser.parse_args()

phase_blocks = {} # map from PS to list of rows that correspond to the csv

compare_df = pd.read_csv(args.csv)

unphased = 0
phase_block_lengths = {} # map from PS to list of [start, stop] of contiguous phaseblocks
last_phase_block = "-10"
for index in len(compare_df):
    phase_set = compare_df["PS"][index]
    if phase_set == "-1":
        continue
    
    pb = phase_blocks.setdefault(phase_set, [])
    pb.append(index)
    pb_lengths = phase_block_lengths.setdefault(phase_set, [])
    pos = compare_df["POS"][index]
    if phase_set == last_phase_block:
        pb_lengths[-1][1] = pos
    else:
        pb_lengths.append([pos, pos+1])

    last_phase_block = phase_block   

for metablock in phase_block_lengths.keys():
    total_span = 0
    min_pos = sys.maxint
    max_pos = 0 
    total_bases = 0
    contiguous_phaseblocks = 0
    for start_stop in phase_block_lengths[metablock]:
        start = start_stop[0]
        stop = start_stop[1]
        total_bases += stop-start
        min_pos = min(min_pos, start)
        max_pos = max(max_pos, stop)
        contiguous_phaseblocks += 1
    total_span = max_pos - min_pos
    print("metablock ",metablock," span ",total_span," total bases ",total_bases," made up of ",contiguous_phaseblocks," contiguous blocks")

