# Investigate structural variation

Requirements:     
- SRA toolkit    


### 01. Obtain raw data files ###
The files are large, and so it is required that `sra-toolkit` is installed. To install, follow the installation instructions on the following [GitHub page](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).          
Another useful resource for [prefetch and fasterq-dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)         

1. Prefetch     
Go to the directory in which you would like to download the files, then run the following command:     
`prefetch SRR23731556`        

This will download a binary containing the required data.    

2. fasterq-dump
`fasterq-dump SRR23731556`    
Be ready though, because this will require a lot of space, for example, a 10 Gb SRA accession (binary) may expand to produce two fastq files of 20 Gb each prior to compression. The website warns that results may be even larger.     

3. Compress
To reduce the size a bit, use gzip:    
`gzip *.fastq` 

