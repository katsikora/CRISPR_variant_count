import os
import sys
import time
import argparse
import re
import logging
import subprocess
import fnmatch

sys.path.insert(0, "/data/boehm/group/pipelines/ruffus")
from ruffus import *
from ruffus.proxy_logger import *
from ruffus.combinatorics import *
from ruffus.drmaa_wrapper import run_job, run_job_using_drmaa, error_drmaa_job

from PIL import Image
import string

parser = argparse.ArgumentParser(prog="CRISPRpipe", version="0.0.1", description="produces sequence clusters from amplicon seq", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ri", "--readIn", dest="readdir", action="store", default=False, help="input read folder")
parser.add_argument("-g", "--ref", dest="refID", action="store", default='GRCm38', help="genome version identifier")
parser.add_argument("-tgt", "--targetInterval", dest="targetInterval", action="store", default=None, help="amplicon genomic interval")
parser.add_argument("-si", "--sampleInfo", dest="sampleInfo", action="store", default=None, help="sample sheet")
parser.add_argument("-wd", "--wdir", dest="wdir", action="store", default=False, help="output folder")
parser.add_argument("-bs", "--batchSize", dest="bsize", action="store",metavar="INT",type=int,default=10, help="number of samples to process in parallel")
parser.add_argument("-nt", "--numThr", dest="nthreads", action="store",metavar="INT",type=int,default=8, help="number of threads to use per sample")
parser.add_argument("-ec", "--errorCorrect", dest="errorCorrect", action="store_true", help="run error correction of sequences (experimental)")
parser.add_argument("--touchOnly", dest="touchOnly", action="store_true", help="only touch files")
parser.add_argument("--target_tasks", dest="target_tasks", action="store",default=[], help="target tasks")
parser.add_argument("--forcedtorun_tasks", dest="forcedtorun_tasks", action="store",default=[], help="forced to run tasks")

args = parser.parse_args()

#setup central working directory
wdir=args.wdir
if not os.path.exists(wdir):
    os.makedirs(wdir)

#setup logging
logger = logging.getLogger(__name__)
fhandler = logging.FileHandler(filename=os.path.join(wdir,'pipeline.log'), mode='a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.debug(subprocess.check_output('echo $DRMAA_LIBRARY_PATH',shell=True))

import drmaa

#initiate 1 drmaa session for the whole pipeline
mySession=drmaa.Session()
mySession.initialize()

#identify pipeline input files
readdir=args.readdir
os.chdir(readdir)

libset2 = []
# Walk through directory
for dName, sdName, fList in os.walk(readdir):
    for fileName in fList:
        if fnmatch.fnmatch(fileName, "*fastq.gz"): # Match search string
            libset2.append(os.path.join(dName, fileName))

libset2_R1=filter(lambda x:'_R1.fastq.gz' in x, libset2)
libset2_R1.sort()
libset2_R2=filter(lambda x:'_R2.fastq.gz' in x, libset2)
libset2_R2.sort()
read_root=[ re.sub('_R1.fastq.gz','',R1f) for R1f in libset2_R1 ]
INfiles=list(zip(libset2_R1,libset2_R2))
    
logger.debug(INfiles)    

######################## paths to unconverted reference genome #######################################
if os.path.isfile(os.path.join('/data/repository/organisms',(args.refID + '_ensembl'),'BWAindex/genome.fa.bwt')):
    refG=os.path.join('/data/repository/organisms',(args.refID + '_ensembl'),'BWAindex/genome.fa')
    
    logger.info('Reference genome ' + refG + 'will be used.' )
else:
    logger.error('Reference genome not recognized.Please check spelling and/or path.')

##################PATHS TO EXECUTABLES###############################################################
FQCpath='/package/FastQC-0.11.3'
cutpath='/package/cutadapt-1.9.1/bin'
prinpath='/data/boehm/sikora/tools/prinseq-lite-0.20.4' #/data/boehm/group/pipelines/tools/COMPLETE_PATH
flashpath='/data/boehm/sikora/tools/FLASH-1.2.11'
bwapath='/package/bwa-0.7.15/bin'
sampath='/package/samtools-1.3/bin'
BEDpath='/package/bedtools2-2.25.0/bin'
Rpath='/package/R-3.3.1/bin'
BBtools='/data/boehm/sikora/tools/bbmap'
FASTXpath='/package/fastx_toolkit_0.0.13'
stkpath='/data/boehm/sikora/tools/seqtk-master'


########adapter-trim and BQ-trim sequencing reads########
cutout=os.path.join(wdir,'reads_cut')
fqcout=os.path.join(wdir,'fastqc_cut')
@mkdir(cutout,os.path.join(cutout,'logs'))
@transform(INfiles,suffix('_R1.fastq.gz'),['_R1.fastq.gz','_R2.fastq.gz'],output_dir=cutout)
def cut_reads(infiles, outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(cutpath,'cutadapt') + ' -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5  -o ' + oo1 + ' -p ' + oo2 + ' ' + ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(cutout,"logs","%s.cut_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.cut_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = bshcmd,
                                      job_name          = 'cut_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("Cut_reads error: %s" % err)
            raise
        else:
           logger.info('Adapter trimming complete')
           
@transform(cut_reads,suffix('_R1.fastq.gz'),['_prin_1.fastq.gz','_prin_2.fastq.gz'],output_dir=cutout)
def trim_BQ(infiles,outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    uzcmd1='zcat -v '+ ii1 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1)))
    uzcmd2='zcat -v '+ ii2 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2)))
    bshcmd='perl '+ os.path.join(prinpath,'prinseq-lite.pl') + ' -fastq ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' -fastq2 ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' -out_good ' + os.path.join(cutout,re.sub('_1.fastq.gz','',os.path.basename(oo1))) +' -trim_qual_right 24 -trim_qual_type mean -trim_qual_window 6 -trim_qual_step 3 -min_len 50 -ns_max_p 10 -out_bad null'
    zcmd1='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' > ' + oo1
    zcmd2='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2))) + ' > ' + oo2
    clcmd='rm -v '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2)))
    cmd_all=';'.join([uzcmd1,uzcmd2,bshcmd,zcmd1,zcmd2,clcmd])
    logger.info(cmd_all)           
    with open(os.path.join(cutout,"logs","%s.BQtrim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.BQtrim_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = cmd_all,
                                      job_name          = 'BQtrim_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("BQtrim_reads error: %s" % err)
            raise
        else:
           logger.info('Base quality trimming complete')
           
@mkdir(fqcout,os.path.join(fqcout,'logs'))            
@transform(trim_BQ,suffix('_1.fastq.gz'),output=['_1.zip','_1.html','_2.zip','_2.html'],output_dir=fqcout)
def postTrim_fqc(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(FQCpath,'fastqc ')+' --outdir ' + fqcout + ' -t 8 '+ ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(fqcout,"logs","%s.post_fqc.out" % read_root),'w+') as stdoutF, open(os.path.join(fqcout,"logs","%s.post_fqc.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'post_fqc',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Post_trim_fastqc error: %s" % err)
            raise
        else:
            logger.info('Post trim fastqc complete')   

#######merge forward and reverse mates########
@transform(trim_BQ,suffix('_1.fastq.gz'),output='_flash.extendedFrags.fastq.gz',output_dir=cutout)
def merge_mates(input_files,output_file):
    ii1 = input_files[0]
    ii2 = input_files[1]
    oo = re.sub('.extendedFrags.fastq.gz','',os.path.basename(output_file))
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(flashpath,'flash')+ ' -z -M 300 -t 8 -o '+ oo + ' -d ' + cutout + ' ' + ii1 + ' ' + ii2 
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.flash.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.flash.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'flash',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Flash error: %s" % err)
            raise
        else:
            logger.info('Merging mates complete')

#### map merged reads to reference genome #####
bamout=os.path.join(wdir,'bams')
@mkdir(bamout,os.path.join(bamout,'logs'))
@transform(merge_mates,suffix('_flash.extendedFrags.fastq.gz'),['.sorted.bam','.sorted.bam.bai'],output_dir=bamout)
def map_reads(input_file,output_files):
    ii=input_file
    oo1=output_files[0]
    read_root=re.sub('_prin_flash.extendedFrags.fastq.gz','',os.path.basename(ii))
    mapcmd=os.path.join(bwapath,'bwa')+' mem -t 8 '+ refG + ' '+ ii + ' | ' + os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -m 6G -@ 8 -o ' + oo1
    indcmd=os.path.join(sampath,'samtools') + ' index ' + oo1
    cmd_all=';'.join([mapcmd,indcmd])
    logger.info(cmd_all)
    with open(os.path.join(bamout,"logs","%s.bwa.out" % read_root),'w+') as stdoutF, open(os.path.join(bamout,"logs","%s.bwa.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd_all,
                                          job_name          = 'bwa',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bwa error: %s" % err)
            raise
        else:
            logger.info('Read mapping complete')

#### output alignments to bed file ################## 
@transform(map_reads,suffix('.sorted.bam'),'.sorted.bed',output_dir=bamout)
def bam2bed(input_files,output_file):
    ii1=input_files[0]
    oo=output_file
    read_root=re.sub('.sorted.bam','',os.path.basename(ii1))
    bshcmd=os.path.join(BEDpath,'bedtools')+' bamtobed -i ' + ii1 + ' > ' + oo
    logger.info(bshcmd)
    with open(os.path.join(bamout,"logs","%s.bam2bed.out" % read_root),'w+') as stdoutF, open(os.path.join(bamout,"logs","%s.bam2bed.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'bam2bed',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bamtobed error: %s" % err)
            raise
        else:
            logger.info('Bam to bed conversion complete')
            
#### filter alignments according to expected amplicon size ####
@transform(bam2bed,suffix('.sorted.bed'),'.whitelist.txt',output_dir=bamout)
def filter_bed(input_file,output_file):
    ii=input_file
    oo=output_file
    read_root=re.sub('.sorted.bed','',os.path.basename(ii))
    Rcmd=os.path.join(Rpath,'Rscript')+' --no-save --no-restore /data/boehm/group/pipelines/CRISPR_variant_counting/v0.0.1/filt_bed.R ' + bamout + ' ' + ii + ' ' + args.targetInterval + ' ;sleep 300'
    logger.info(Rcmd)
    with open(os.path.join(bamout,"logs","%s.bed_filt.out" % read_root),'w+') as stdoutF, open(os.path.join(bamout,"logs","%s.bed_filt.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = Rcmd,
                                          job_name          = 'bed_filt',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bed filtering error: %s" % err)
            raise
        else:
            logger.info('Bed filtering complete')
            
#### filter bam file for whitelist read names; output to fasta #####################
seqout=os.path.join(wdir,'filtered_sequences')
@follows(mkdir(seqout),mkdir(os.path.join(seqout,"logs")))
@transform(input=filter_bed,filter=suffix('.whitelist.txt'),output='.whitelist.fastq.gz',output_dir=seqout)
def filter_bam(input_file,output_file):
    ii1=input_file
    ii2=re.sub('.whitelist.txt','.sorted.bam',ii1)
    oo=output_file
    read_root=re.sub('.whitelist.txt','',os.path.basename(ii1))
    cmd=os.path.join(BBtools,'filterbyname.sh')+ ' -Xmx20g in=' + ii2 + ' out=' + oo + ' names='+ ii1 + ' include=true overwrite=true; sleep 300'
    logger.info(cmd)
    with open(os.path.join(seqout,"logs","%s.bam_filt.out" % read_root),'w+') as stdoutF, open(os.path.join(seqout,"logs","%s.bam_filt.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd,
                                          job_name          = 'bam_filt',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=20000')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bam filtering error: %s" % err)
            raise
        else:
            logger.info('Bam filtering complete')

@transform(input=filter_bam,filter=suffix('.whitelist.fastq.gz'),output='.whitelist.fasta',output_dir=seqout)
def fq2fa(input_file,output_file):
    ii=input_file
    oo=output_file
    read_root=re.sub('.whitelist.fastq.gz','',os.path.basename(ii))
    cmd='zcat -v '+ ii + ' | awk \'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}\' - > ' + oo + ' ;sleep 300'
    logger.info(cmd)
    with open(os.path.join(seqout,"logs","%s.fq2fa.out" % read_root),'w+') as stdoutF, open(os.path.join(seqout,"logs","%s.fq2fa.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd,
                                          job_name          = 'fq2fa',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Fastq to fasta conversion error: %s" % err)
            raise
        else:
            logger.info('Fastq to fasta conversion complete')
         
            
#### cluster identical sequences of the same length in the fasta file ######
@transform(fq2fa,suffix('.whitelist.fasta'),'.whitelist.dedup.fasta',output_dir=seqout)
def fa_clust(input_file,output_file):
    ii=input_file
    oo=output_file
    read_root=re.sub('.whitelist.fasta','',os.path.basename(ii))
    bshcmd=os.path.join(FASTXpath,'fastx_collapser')+ ' -i ' +ii + ' -o ' + oo + ' ; sleep 300'
    logger.info(bshcmd)
    with open(os.path.join(seqout,"logs","%s.fa_clust.out" % read_root),'w+') as stdoutF, open(os.path.join(seqout,"logs","%s.fa_clust.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'fa_clust',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Fasta clustering error: %s" % err)
            raise
        else:
            logger.info('Fasta clustering complete')

################# postprocess cluster counts in R #########################################
Rout=os.path.join(wdir,'workspaceR')
if not os.path.exists(Rout):
    os.makedirs(Rout)
os.chdir(Rout)
@follows(mkdir(Rout),mkdir(os.path.join(Rout,'logs')))
@merge(fa_clust,output=os.path.join(Rout,'Cluster.counts.txt'))
def R_get_cluster_counts(input_files,output_file):
    ii=input_files[1]
    idir=seqout
    oo=output_file
    read_root=re.sub('.whitelist.dedup.fasta','',os.path.basename(ii))
    Rcmd=os.path.join(Rpath,'Rscript')+' --no-save --no-restore /data/boehm/group/pipelines/CRISPR_variant_counting/v0.0.1/cluster_counts.R ' + Rout + ' ' + idir
    logger.info(Rcmd)
    with open(os.path.join(Rout,"logs","%s.clu_counts.out" % read_root),'w+') as stdoutF, open(os.path.join(Rout,"logs","%s.clu_counts.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = Rcmd,
                                          job_name          = 'clu_counts',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Cluster counting error: %s" % err)
            raise
        else:
            logger.info('Cluster counting complete')
         
         
#####main

if __name__ == '__main__':
    with open(os.path.join(wdir,"pipelineGraph.png"),'w') as pipeGraph:
        pipeline_printout_graph(stream=pipeGraph,output_format='png',pipeline_name='WGBS',target_tasks=args.target_tasks)
    with open (os.path.join(wdir,"pipelinePrint.txt"),'w') as pipePrint:
        pipeline_printout(verbose_abbreviated_path=0,output_stream=pipePrint,target_tasks=args.target_tasks)    

    pipeline_run(touch_files_only=args.touchOnly,multiprocess=args.bsize,target_tasks=args.target_tasks,forcedtorun_tasks=args.forcedtorun_tasks,logger=logger)
    mySession.exit()         


