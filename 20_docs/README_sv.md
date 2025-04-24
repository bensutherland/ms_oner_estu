# Investigate structural variation

Requirements:     
- SRA toolkit    
- E. Normandeau 'scripts' repository
- github repository wgrs_workflow

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

This will operate on the folder containing the .sra file, and will put the fastq files directly into the directory from which the command was made (generally above the directory holding the .sra).   

3. Compress
To reduce the size a bit, use gzip:    
`gzip *.fastq` 


### 02. Obtain and prepare reference genome ###
Download the reference genome from NCBI [GCF_034236695.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034236695.1/).       

Index with bwa     
`bwa index GCF_034236695.1_Oner_Uvic_2.0_genomic.fna`      

Extract only chromosome 18       
```
# Copy genome to repository
cp -l ~/genomes/GCF_034236695.1_Oner_Uvic_2.0_genomic.fna ./04_genome/ 

# Change into genome folder
cd 04_genome

# Unwrap genome
fasta_unwrap.py ./GCF_034236695.1_Oner_Uvic_2.0_genomic.fna GCF_034236695.1_Oner_Uvic_2.0_genomic_unwrap.fna

# Obtain chr18 only
grep -A1 'NC_088413' GCF_034236695.1_Oner_Uvic_2.0_genomic_unwrap.fna > GCF_034236695.1_Oner_Uvic_2.0_genomic_unwrap_chr18_only.fna

# index chr18 genome
bwa index GCF_034236695.1_Oner_Uvic_2.0_genomic_unwrap_chr18_only.fna
```

### 03. Trim and align raw data to subset reference genome ###
Clone the github repository `wgrs_workflow`. Follow README.       


