$ cd /work/Neuroinformatics_Core/s204365/ATACSeq_0001/metadata
$ module load dos2unix/gcc/7.3.4

# Data inconsistency (if utilizing Windows OS) caused by preparing the data in Windows and passing it onto Linux. There is often an "M" character that is falsely recognized in Linux. 
# Remove this special character from data file(s) by:
$ dos2unix fragments_104.tsv.gz fragments_105.tsv.gz fragments_106.tsv.gz fragments_107.tsv.gz

$ cut -d " " -f 2- fragments_104.tsv.gz > fragments104_processed.tsv.gz
$ cut -d " " -f 2- fragments_105.tsv.gz > fragments105_processed.tsv.gz
$ cut -d " " -f 2- fragments_106.tsv.gz > fragments106_processed.tsv.gz
$ cut -d " " -f 2- fragments_107.tsv.gz > fragments107_processed.tsv.gz

