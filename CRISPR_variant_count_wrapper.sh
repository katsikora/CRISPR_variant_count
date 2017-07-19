DRMAA_LIBRARY_PATH=/data/manke/repository/scripts/DNA_methylation/drmaa/lib/libdrmaa.so
export PATH=$PATH:/package/bwa-0.7.4/bin
export R_LIBS_USER=/data/manke/repository/scripts/DNA_methylation/Rlibs.3.3.1
source /data/boehm/sikora/miniconda3/bin/activate NGSpy2.7



python CRISPR_variant_count.py --readIn /data/boehm/sequencing_data/170124_M01358_0018_000000000-B2JYG/Project_A737_Boehm_Diekhoff --ref GRCz10 --targetInterval X:53001978-53002368 --wdir /data/processing3/sikora/iwanami/A737.test --touchOnly

source deactivate

echo 'done all'
