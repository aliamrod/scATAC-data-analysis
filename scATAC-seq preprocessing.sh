$ cd /work/Neuroinformatics_Core/s204365/ATACSeq_0001/metadata
$ module load dos2unix/gcc/7.3.4
$ dos2unix fragments_104.tsv.gz

$ cut -d " " -f 2- fragments_104.tsv.gz > fragments104_processed.tsv.gz
$ cut -d " " -f 2- fragments_105.tsv.gz > fragments105_processed.tsv.gz
$ cut -d " " -f 2- fragments_106.tsv.gz > fragments106_processed.tsv.gz
$ cut -d " " -f 2- fragments_107.tsv.gz > fragments107_processed.tsv.gz

# 
