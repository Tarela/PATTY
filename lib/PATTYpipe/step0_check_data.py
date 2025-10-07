#!/usr/bin/env python
 
# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *
import platform

# --------------------------
# custom package
# --------------------------

### tool function

### tool function
import PATTYpipe
from PATTYpipe.Utility import (sp, 
                    raise_error,
                    wlog,
                    ewlog,
                    checkbedformat,
                    CMD
                    )

# --------------------------
# main 
# --------------------------

def step0_check_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    ### check data path , format ,
    # input cuttag&atac dat
    if "/" in conf_dict['General']['cuttag']:
        if conf_dict['General']['cuttag'].startswith("/"):
            pass
        elif conf_dict['General']['cuttag'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['cuttag'] = homedir +"/" + conf_dict['General']['cuttag'][1:]
        else:
            conf_dict['General']['cuttag'] = conf_dict['General']['startdir'] + conf_dict['General']['cuttag']
    else:
        conf_dict['General']['cuttag'] = conf_dict['General']['startdir'] + conf_dict['General']['cuttag']
    if not os.path.isfile(conf_dict['General']['cuttag']):
        ewlog("cuttag file %s not found"%(conf_dict['General']['cuttag']),logfile)
    if not conf_dict['General']['cuttag'].endswith('.bed') and not conf_dict['General']['cuttag'].endswith('.bed.gz'):
        ewlog('extenion of cuttag file is not .bed (nor .bed.gz)',logfile)
    checkbed = checkbedformat(conf_dict['General']['cuttag'])
    if checkbed == "fail":
        ewlog("cuttag file is not a PE/SE bed file",logfile)
    else:#elif checkbed in ["PE","SE"]:
        conf_dict['General']['format'] = checkbed
        wlog("detected bed file format is %s"%checkbed,logfile)

    if "/" in conf_dict['General']['atac']:
        if conf_dict['General']['atac'].startswith("/"):
            pass
        elif conf_dict['General']['atac'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['atac'] = homedir +"/" + conf_dict['General']['atac'][1:]
        else:
            conf_dict['General']['atac'] = conf_dict['General']['startdir'] + conf_dict['General']['atac']
    else:
        conf_dict['General']['atac'] = conf_dict['General']['startdir'] + conf_dict['General']['atac']
    if not os.path.isfile(conf_dict['General']['atac']):
        ewlog("atac file %s not found"%(conf_dict['General']['atac']),logfile)
    if not conf_dict['General']['atac'].endswith('.bed') and not conf_dict['General']['atac'].endswith('.bed.gz'):
        ewlog('extenion of atac file is not .bed (nor .bed.gz)',logfile)
    checkbed = checkbedformat(conf_dict['General']['atac'])
    if checkbed == "fail":
        ewlog("atac file is not a PE/SE bed file",logfile)
    else:#elif checkbed in ["PE","SE"]:
        wlog("detected bed file format is %s"%checkbed,logfile)

    if "/" in conf_dict['General']['genomebin']:
        if conf_dict['General']['genomebin'].startswith("/"):
            pass
        elif conf_dict['General']['genomebin'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['genomebin'] = homedir +"/" + conf_dict['General']['genomebin'][1:]
        else:
            conf_dict['General']['genomebin'] = conf_dict['General']['startdir'] + conf_dict['General']['genomebin']
    else:
        conf_dict['General']['genomebin'] = conf_dict['General']['startdir'] + conf_dict['General']['genomebin']
    if not os.path.isfile(conf_dict['General']['genomebin']):
        ewlog("genome-wide mappable bin file %s not found"%(conf_dict['General']['genomebin']),logfile)
    if not conf_dict['General']['atac'].endswith('.bed') and not conf_dict['General']['genomebin'].endswith('.bed.gz'):
        ewlog('extenion of genome-wide mappable bin file is not .bed (nor .bed.gz)',logfile)

    # --binlist
    if conf_dict['options']['binlist']:
        if "/" in conf_dict['options']['binlist']:
            if conf_dict['options']['binlist'].startswith("/"):
                pass
            elif conf_dict['options']['binlist'].startswith("~/"):
                homedir = os.path.expanduser("~")
                conf_dict['options']['binlist'] = homedir +"/" + conf_dict['options']['binlist'][1:]
            else:
                conf_dict['options']['binlist'] = conf_dict['General']['startdir'] + conf_dict['options']['binlist']
        else:
            conf_dict['options']['binlist'] = conf_dict['General']['startdir'] + conf_dict['options']['binlist']

        if not os.path.isfile(conf_dict['options']['binlist']):
            wlog("external binlist file %s not found, PATTY use genome-wide mappable bins with sufficient reads"%(conf_dict['options']['binlist']),logfile)
            conf_dict['options']['binlist'] = "NA"
        checkbed = checkbedformat(conf_dict['options']['binlist'])
        if checkbed == "fail":
            wlog("external binlist file %s is not a bed file, PATTY use genome-wide mappable bins with sufficient reads"%(conf_dict['options']['binlist']),logfile)
            conf_dict['options']['binlist'] = "NA"
    else:
        wlog("no external binlist file inputted, PATTY use genome-wide mappable bins file as target regions",logfile)
        conf_dict['options']['binlist'] = "NA"
        
    # check software
    OS = platform.system()
    # check system
    if OS == "Linux":
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"
    elif OS == "Darwin":
        bedtools_software = "bedtools_mac"
        bdg2bw_software = "bedGraphToBigWig_mac"
    else:
        wlog("detected system is nither linux nor mac, try linux version",logfile)
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"

    check_bedtools = sp("which bedtools")
    check_bdg2bw = sp("which bedGraphToBigWig")
    check_R = sp("which Rscript")
#    check_pdflatex = sp("which pdflatex")

    if check_bedtools[0].decode("ascii") != "":
        wlog("bedtools installed",logfile)
        conf_dict['General']['bedtools'] = "bedtools"
    else:
        wlog("bedtools not installed in the default path, use PATTY-built-in bedtools",logfile)
        conf_dict['General']['bedtools'] = PATTYpipe.__path__[0]+"/external_script/%s"%bedtools_software

    if check_bdg2bw[0].decode("ascii") != "":
        wlog("bedGraphToBigWig(UCSCtools) installed",logfile)
        conf_dict['General']['bedGraphToBigWig'] = "bedGraphToBigWig"
    else:
        wlog("bedGraphToBigWig(UCSCtools) not installed in the default path, use PATTY-built-in bedGraphToBigWig",logfile)
        conf_dict['General']['bedGraphToBigWig'] = PATTYpipe.__path__[0]+"/external_script/%s"%bdg2bw_software

    if check_R[0].decode("ascii") != "":
        wlog("Rscript installed",logfile)
        conf_dict['General']['Rscript'] = "Rscript"
    else:
        ewlog("require Rscript installed in the default path",logfile)

    if conf_dict['General']['mode'] == "sc":
#        if conf_dict['options']['topDim'] < 30:
#            wlog("topDim in single-cell clustering should be >=30, set topDim=30",logfile)
#            conf_dict['options']['topDim']=30

        # --cellnames
        if conf_dict['options']["cellnames"]:
            if "/" in conf_dict['options']["cellnames"]:
                if conf_dict['options']["cellnames"].startswith("/"):
                    pass
                elif conf_dict['options']["cellnames"].startswith("~/"):
                    homedir = os.path.expanduser("~")
                    conf_dict['options']["cellnames"] = homedir +"/" + conf_dict['options']["cellnames"][1:]
                else:
                    conf_dict['options']["cellnames"] = conf_dict['General']['startdir'] + conf_dict['options']["cellnames"]
            else:
                conf_dict['options']["cellnames"] = conf_dict['General']['startdir'] + conf_dict['options']["cellnames"]
            if os.path.isfile(conf_dict['options']["cellnames"]):
                conf_dict['options']["usecells"] = []
                wlog("readin used cellnames",logfile)
                inf = open(conf_dict['options']["cellnames"])
                for line in inf:
                    ll = line.strip().split("\t")
                    conf_dict['options']["usecells"].append(ll[0])
                inf.close()
                if len(conf_dict['options']["usecells"]) < 100:
                    wlog("less than 100 cells specified in cellnames file (%s cells). ignore --cellnames parameter"%(len(conf_dict['options']["usecells"])),logfile)
                    conf_dict['options']["usecells"] = []
            else:
                wlog("cellnames file %s not found, ignore --cellnames parameter"%(conf_dict['options']["cellnames"]),logfile)
                conf_dict['options']["usecells"] = []

    ### check trained model file
    conf_dict['General']['model'] = PATTYpipe.__path__[0]+"/refdata/%s_LR_CnTATAC_model.joblib"%(conf_dict['General']["factor"])

    return conf_dict

#    conf_dict['options']['chromosome'] = []
#    inf = open(conf_dict['options']['csize'])
#    for line in inf:
#        if not line.strip()=="":
#            conf_dict['options']['chromosome'].append(line.strip().split("\t")[0])
#    inf.close()

    ### chromosome list
#    if conf_dict['General']["genome"] == "hg38":
#        conf_dict['General']["chrmlist"] = []
#    else:
#        conf_dict['General']["chrmlist"] = []

    #if not os.path.isfile(conf_dict['options']['biasfile']):
    #    wlog("no naked DNA bias matrix, use mtDNA reads to estimate",logfile)
    #    conf_dict['options']['bias'] = "chrM"


    ### check pdflatex
    #if sp('pdflatex --help')[0] == "":
    #    wlog('pdflatex was not installed, ncHMR_detector is still processing but no summary report generated',logfile)
    #    conf_dict['General']['latex'] = 0
    #else:
    #    conf_dict['General']['latex'] = 1


 #   if conf_dict['options']['signalname']

#    conf_dict['General']['signalname'] = []
#    conf_dict['General']['signalfile'] = []
#    for bwsignalfile in conf_dict['General']['signal']:
#        if not bwsignalfile.startswith('/'):
#            bwsignalfile = conf_dict['General']['startdir'] + bwsignalfile
#
#        if not os.path.isfile(bwsignalfile):
#            wlog("signal bw file %s not found, ignored"%(bwsignalfile),logfile)
#            continue
#
#        if bwsignalfile.endswith(".bw"):
#            conf_dict['General']['signalfile'].append(bwsignalfile)
#            conf_dict['General']['signalname'].append(bwsignalfile.split("/")[-1][:-3])
#        elif bwsignalfile.endswith(".bigwig"):
#            conf_dict['General']['signalfile'].append(bwsignalfile)
#            conf_dict['General']['signalname'].append(bwsignalfile.split("/")[-1][:-7])
#        else:
#            wlog('[WARNING] extension of signal bw file is not bw/bigwig',logfile)
#            conf_dict['General']['signalfile'].append(bwsignalfile)
#            conf_dict['General']['signalname'].append(bwsignalfile.split("/")[-1])
#
#    if len(conf_dict['General']['signalfile']) == 0:
#        ewlog("no signal bw file valid, exit")
#    elif len(conf_dict['General']['signalfile']) > 4:
#        ewlog("maximum signal bw file is limited to 4. There were %s signal file inputed, exit"%(len(conf_dict['General']['signalfile'])))
#
#    ### check TFpeak folder    
#    if "~" in conf_dict['General']['peakFolder']:
#        ewlog('require absolute path for peak/track Folder, Folder cannot contain "~", current Folder is %s'%(conf_dict['General']['peakFolder']),logfile)
#    if not conf_dict['General']['peakFolder'].startswith('/'):
#        conf_dict['General']['peakFolder'] = conf_dict['General']['startdir'] + conf_dict['General']['peakFolder']
#    if not conf_dict['General']['peakFolder'].endswith('/'):
#        conf_dict['General']['peakFolder'] += "/"
#    if not os.path.isdir(conf_dict['General']['peakFolder']):
#        ewlog("Folder %s not found"%(conf_dict['General']['peakFolder']),logfile)
#
#    if conf_dict['General']['mode'] == "signal":
#        wlog("signal mode is activated",logfile)
#        if conf_dict['General']['bwfolder']:
#            wlog("bwFolder is specified, checking data for signal mode",logfile)
#            if "~" in conf_dict['General']['bwfolder']:
#                wlog('require absolute path for bwFolder, bwFolder cannot contain "~", current Folder is %s, use peak mode'%(conf_dict['General']['bwfolder']),logfile)
#                conf_dict['General']['mode'] = "binary"
#            else:
#                if not conf_dict['General']['bwfolder'].startswith('/'):
#                    conf_dict['General']['bwfolder'] = conf_dict['General']['startdir'] + conf_dict['General']['bwfolder']
#                if not conf_dict['General']['bwfolder'].endswith('/'):
#                    conf_dict['General']['bwfolder'] += "/"
#                if not os.path.isdir(conf_dict['General']['bwfolder']):
#                    wlog("bwFolder %s not found, use binary mode"%(conf_dict['General']['peakFolder']),logfile)
#                    conf_dict['General']['mode'] = "binary"
#        else:
#            wlog("bwfolder is not specified, use binary mode",logfile)
#            conf_dict['General']['mode'] = "binary"
#    else:
#        wlog("binary mode is activaed",logfile)
#
#
#    wlog("Check the peak.bed files in the Folder, only '.bed' files with >1000 peaks are included in the following analysis",logfile)
#    conf_dict['General']['peakfilenames'] = []
#    for f in os.listdir(conf_dict['General']['peakFolder']):
#        if f.endswith(".bed") and os.path.isfile(conf_dict['General']['peakFolder']+f):
#            checkbed = checkbedformat(conf_dict['General']['cuttag'],1000)
#            if checkbed == "pass":
#                conf_dict['General']['peakfilenames'].append(f[:-4])
#
#    if (len(conf_dict['General']['peakfilenames']) == 0):
#        ewlog("no peak file (cofactor candidate) in (bed format & >1000peaks) are included, exit",logfile)
#
#    if conf_dict['General']['mode'] == "signal":
#        conf_dict['General']['bwfilenames'] = []
#        for f in os.listdir(conf_dict['General']['bwfolder']):
#            if f.endswith(".bw") and os.path.isfile(conf_dict['General']['bwfolder']+f):
#                conf_dict['General']['bwfilenames'].append(f[:-3])
#        ### compare the name from bwfiles and peak files
#        conf_dict['General']['usefilename'] = []
#        for name in conf_dict['General']['bwfilenames']:
#            if name in conf_dict['General']['peakfilenames']:
#                conf_dict['General']['usefilename'].append(name)
#        ### if less than 50% peakfiles share name with bwfiles, change back to peak mode
#        if len(conf_dict['General']['usefilename']) < len(conf_dict['General']['peakfilenames'])*0.5:
#            conf_dict['General']['mode'] = "binary"
#            wlog("the number of shared peak&bw files is less than half of the number of peakfiles, use binary mode",logfile)
#        else:
#            wlog("all checks for signal mode passed, use signal mode",logfile)
#
#    if conf_dict['General']['mode'] == "binary":
#        conf_dict['General']['usefilename'] = conf_dict['General']['peakfilenames']
#
#    wlog("%s cofactor candidates are included"%(len(conf_dict['General']['usefilename'])),logfile)
#
#                #checkbed = checkbedformat(conf_dict['General']['cuttag'],1000)
#                #if checkbed == "pass":
#                #conf_dict['General']['peakfilenames'].append(f[:-3])            
#
#    outf = open(conf_dict['General']['outname']+"_cofactor_candidate_list.txt",'w')
#    for cofactor in conf_dict['General']['usefilename']:
#        outf.write(cofactor+"\n")
#    outf.close()
#    ### check options
#    wlog('check option: ',logfile)
##
#    try:
#        wlog("extend length for HMsignal is %s bp"%(int(conf_dict['options']['ext'])),logfile)
#        conf_dict['options']['ext']=int(conf_dict['options']['ext'])
#    except:
#        wlog("extend length %s is not valid, use default value: 1000bp"%(conf_dict['options']['ext']),logfile)
#        conf_dict['options']['ext'] = 1000
#
#    try:
#        wlog("use Pvalue = %s as cutoff"%(str(float(conf_dict['options']['Pvalue']))),logfile)
#    except:
#        wlog("input Pvalue %s is not recognized, use default Pvalue=0.001"%(conf_dict['options']['Pvalue']),logfile)
#        conf_dict['options']['Pvalue'] = 0.001 
#    if float(conf_dict['options']['Pvalue']) >= 1:
#        wlog("input Pvalue %s is not valid, use default Pvalue=0.001"%(conf_dict['options']['Pvalue']),logfile)
#        conf_dict['options']['Pvalue'] = 0.001 
#
#    try: 
#        usealpha = float(conf_dict['options']['Alpha'])
#        if usealpha >=1:
#            wlog("alpha (for elastic-net) cannot be >=1, use alpha=0.5",logfile)
#            conf_dict['options']['Alpha'] = 0.5
#        else:
#            wlog("Alpha (for elastic-net) = %s"%(str(float(conf_dict['options']['Alpha']))),logfile)
#            conf_dict['options']['Alpha']=usealpha
#    except:
#        wlog("input alpha (for elastic-net) %s is not valid, use alpha=0.5"%(conf_dict['options']['Alpha']),logfile)
#
#    wlog("Lambda choice is %s"%(conf_dict['options']['Lambda']),logfile)
#    if conf_dict['options']['TopNcofactors'] == "all":
#        wlog("all significant co-factors will be output",logfile)
#    else:
#        try:
#            topTF = int(conf_dict['options']['TopNcofactors'])
#            wlog("the topN number %s will be output"%(conf_dict['options']['TopNcofactors']),logfile)
#            conf_dict['options']['TopNcofactors'] = topTF
#        except:
#            wlog("the topN number %s is not valid, output top5 co-factors"%(conf_dict['options']['TopNcofactors']),logfile)
#            conf_dict['options']['TopNcofactors'] = 5
#

    
    
