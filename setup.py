#!/usr/bin/env python
"""Description
Setup script for "PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data"
Copyright (c) 2025 Shengen Hu <sh8tv@virginia.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
""" 
import os
import sys
import subprocess
import platform
from distutils.core import setup#, Extension
import distutils.command.install_lib
#from setuptools import setup, find_packages

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1   


class my_install_lib(distutils.command.install_lib.install_lib):
    def run(self):
        distutils.command.install_lib.install_lib.run(self)
        mode = 755
        # here we start with doing our overriding and private magic ..
        for filepath in self.get_outputs():
            if  "bedtools" in filepath or "bedGraphToBigWig" in filepath or "twoBitToFa" in filepath or "twoBitInfo" in filepath:
            #if self.install_scripts in filepath:
            #    log.info("Overriding setuptools mode of scripts ...")
            #    log.info("Changing ownership of %s to uid:%s gid %s" %
            #             (filepath, uid, gid))
            #    os.chown(filepath, uid, gid)
            #    log.info("Changing permissions of %s to %s" %
            #             (filepath, oct(mode)))
                os.chmod(filepath, mode)

def main(): 
#    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
#	    print >> sys.stderr, "ERROR: ncHMR_detector requires Python 2.7"
#	    sys.exit()
    has_R = check_R()
    if has_R == 0:
	    print("ERROR: PATTY requires R & Rscript under default PATH", file=sys.stderr)
	    sys.exit()

    OS = platform.system()
#    if OS == "Linux":
#        bwsum_software = "bigWigSummary_linux"
#    elif OS == "Darwin":
#        bwsum_software = "bigWigSummary_mac"
#    else:
#        wlog("detected system is nither linux nor mac, try linux version of bigWigSummary",logfile)
#        bwsum_software = "bigWigSummary_linux"

    OS = platform.system()
    if OS == "Linux":
#        bwsum_software = "bigWigSummary_linux"
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"
#       twobit_software = "twoBitToFa_linux"
#       twobitI_software = "twoBitInfo_linux"
    elif OS == "Darwin":
#        bwsum_software = "bigWigSummary_mac"
        bedtools_software = "bedtools_mac"
        bdg2bw_software = "bedGraphToBigWig_mac"
#       twobit_software = "twoBitToFa_mac"
#       twobitI_software = "twoBitInfo_mac"
    else:
        wlog("detected system is nither linux nor mac, try linux version",logfile)
#        bwsum_software = "bigWigSummary_linux"
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"
#        twobit_software = "twoBitToFa_linux"
#        twobitI_software = "twoBitInfo_linux"
                
    setup(name="PATTY",
          version="1.0",
          description="PATTY: a computational method for correcting open chromatin bias in bulk and single-cell CUT&Tag data",
          author='Shengen Shawn Hu',
          author_email='sh8tv@virginia.edu',
          url='https://github.com/Tarela/PATTY.git',
          package_dir={'PATTYpipe' : 'lib'},
          packages=['PATTYpipe'],
          #package_data={}
          package_data={'PATTYpipe': [#'external_script/%s'%bwsum_software,
                                      'external_script/%s'%bedtools_software,
                                      'external_script/%s'%bdg2bw_software,
                                      #'external_script/%s'%twobit_software,
                                      #'external_script/%s'%twobitI_software,
                                      'refdata/H3K27me3_LR_CnTATAC_model.joblib',
                                      'refdata/H3K27ac_LR_CnTATAC_model.joblib',
                                      'refdata/H3K9me3_LR_CnTATAC_model.joblib',
                                      'refdata/hg38_36_chr_fareadsClean_ge100mapregion.bed.gz',
                                      'refdata/mm10_36_chr_fareadsClean_ge100mapregion.bed.gz'
                                      ]},#,#'Config/template.conf',
                                  #'Rscript/analysis.r',
                                  #'Rscript/individual_qc.r',
                                  #'Rscript/readsbulkQC.r',
                                  #'Rscript/detectNonCanonical.r'
                                  #   ]},
          #scripts=['bin/ncHMR_detector_py3','refpackage/bwsummary/%s'%bwsum_software],
          scripts=['bin/PATTY'],#,'refpackage/bwsummary/%s'%bwsum_software],
          #data_files=[('/Users/sh8tv/bin',['refpackage/bwsummary/%s'%bwsum_software])],
          classifiers=[
        'Development Status :: version1.0 finish',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: pipeline',
        ],
          requires=[],
          cmdclass={'install_lib':my_install_lib}
          )

    print('Installation of PATTY is DONE')


if __name__ == '__main__':
    main()




