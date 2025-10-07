#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
import gzip
import numpy as np
import copy,random
from joblib import load


# --------------------------
# custom package
# --------------------------

### tool function
from PATTYpipe.Utility      import (sp, 
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   open_bed_file,
                                   add_coverPos,
                                   read_in_reads,
                                   prepare_chrom_size
                                   )

# -------------------------- 
# main 
# --------------------------
def step1_BULKscanSig(conf_dict,logfile):

    ### read-in bins
    wlog('read-in bins',logfile)
    # create 200bp bin
    if conf_dict['options']['binlist'] == "NA":
        binFile = conf_dict['General']['genomebin']
    else:
        cmd = """%s intersect -a %s -b %s -u > tmp_usebins.bed"""%(conf_dict['General']['bedtools'],conf_dict['General']['genomebin'],conf_dict['options']['binlist'])
        tmplog = sp(cmd)
        binFile = "tmp_usebins.bed"

    ## readin bins
    binInfo = {}
    binSig = {}
    binextsize=1000
    with open_bed_file(binFile) as inf:
        for line in inf:
            ll = line.split()
            binname = ll[0]+":"+ll[1]+"-"+ll[2]
            binmid = int(ll[1])+ 100
            ext_left = binmid - int(binextsize/2)
            ext_right = binmid + int(binextsize/2)
            if ext_left >= 0:
                binInfo[binname] = ll
                for this_bin in range(ext_left, ext_right, 100):
                    relatedBinname = ll[0]+":"+str(this_bin)+"-"+str(this_bin+100)
                    binSig[relatedBinname] = np.array([0]*100)

    binSig_CnT = copy.deepcopy(binSig)
    binSig_ATAC = copy.deepcopy(binSig)
    wlog("read-in CUT&Tag reads",logfile)
    [binSig_CnT, CnT_total]  = read_in_reads(conf_dict['General']['cuttag'], binSig_CnT)
    wlog("read-in ATAC reads",logfile)
    [binSig_ATAC, ATAC_total] = read_in_reads(conf_dict['General']['atac'],   binSig_ATAC)

    # print bin with siganl
    wlog("bulk correction",logfile)
    model = load(conf_dict['General']['model'])

    outf = open('%s_correctSig.bdg'%(conf_dict['General']['outname']),'w')
    for this_bin in sorted(binInfo.keys()):
        binmid = int(binInfo[this_bin][1])+ 100
        ext_left = binmid - int(binextsize/2)
        ext_right = binmid + int(binextsize/2)
        outValue_CnT = []
        outValue_ATAC = []
        for related_bin in range(ext_left, ext_right, 100):
            relatedBinname = binInfo[this_bin][0]+":"+str(related_bin)+"-"+str(related_bin+100)
            this_array_CnT = binSig_CnT[relatedBinname] / CnT_total * 1e6
            ave_array_CnT = np.mean(this_array_CnT.reshape(10, 10), axis=1)
            this_array_ATAC = binSig_ATAC[relatedBinname] / ATAC_total * 1e6
            ave_array_ATAC = np.mean(this_array_ATAC.reshape(10, 10), axis=1)
            outValue_CnT.extend(list(ave_array_CnT))
            outValue_ATAC.extend(list(ave_array_ATAC))

        if sum(outValue_CnT) == 0:
            continue

        vals_C = np.asarray(outValue_CnT,  dtype=np.float32)  # length 100
        vals_A = np.asarray(outValue_ATAC, dtype=np.float32)  # length 100
        X_row = np.empty((1, 200), dtype=np.float32)
        X_row[0, 0::2] = vals_C
        X_row[0, 1::2] = vals_A
        # Concatenate into a single row (shape: 1 × 200) and predict
        #X_row = np.concatenate([outValue_CnT, outValue_ATAC], axis=0).reshape(1, -1)
        try:
            y_pred = float(model.predict_proba(X_row)[:, 1][0])
        except Exception as e:
            # If your model expects a different shape/order, log and skip
            wlog(f"predict_proba failed for {this_bin}: {e}", logfile)
            continue

        # Write one line: chr  start  end  score
        chr_, start_, end_ = binInfo[this_bin][:3]
        outf.write(f"{chr_}\t{start_}\t{end_}\t{y_pred}\n")

    outf.close()

    wlog("generate bigwig track",logfile)
    prepare_chrom_size(conf_dict['General']['genome'],"tmp_genome.len")
    cmd1 = """sort -k1,1 -k2,2n %s_correctSig.bdg > tmp_sorted.bdg"""%(conf_dict['General']['outname'])
    cmd2 = """%s tmp_sorted.bdg tmp_genome.len %s_correctSig.bw"""%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname'])
    tmplog = sp(cmd1)
    tmplog = sp(cmd2)

    return conf_dict

