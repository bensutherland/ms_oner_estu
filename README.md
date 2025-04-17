# Early Stuart sockeye salmon population genetics
Code to support analysis of Early Stuart sockeye metapopulation within the context of the upper Fraser River sockeye.     

### Whole-genome resequencing analysis ###
#### Requirements ####
- bcftools    
- simple_pop_stats

Commands are run from individual repositories, as indicated in the section.     
Clone `simple_pop_stats`.     

#### Data sources ####
Download source VCF file [Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz](https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Christensen_et_al_2024/25705428) and put it in `simple_pop_stats/02_input_data`.        

Change directory into `simple_pop_stats`.      

#### Subset samples to target region ####
Keep only samples from populations of interest from the BCF file:    
```
# Upper Fraser River samples (no kokanee)
bcftools query -l 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz | grep -E 'Bivouac|Driftwood|Dust|Felix|Paula|Takla|Kuzkwa|Middle|Pinchi|Tachie|Stellako|Nadina|Bowron|Horsefly|BlueLead|McKinley|Mitchell|horsefly|Wasko|Quesnel|Chilko|Taseko' - | grep -vE '\_K\_' - > 02_input_data/samples_to_retain.txt

# Use sample file to subset the VCF file to only these samples 
bcftools view -S 02_input_data/samples_to_retain.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf    

```

#### Filter variants for population genetic use ####
Additional filters to pre-filtered input file:       
```
# Add tags to conduct MAF filter
bcftools +fill-tags 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_w_tags.bcf

# MAF filter
bcftools view -i 'MAF > 0.05' 02_input_data/*_w_tags.bcf -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_w_tags_MAF0.05.bcf

# Generate an LD filtered, compressed VCF file for population genetic analyses     
bcftools +prune -w 50kb -m 0.5 -Oz -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_w_tags_MAF0.05_5w50kb.vcf.gz 02_input_data/*_MAF0.05.bcf

```

#### Population genetic analysis, full region ####
Use the following script to analyze the general population genetic trends in the region:    
`01_scripts/popgen_01_upper_fraser_all.R`       

And to inspect specific regions:    
`01_scripts/popgen_02_upper_fraser_regional.R`    


### Closer inspection of specific region ###
#### All populations, all loci, one chromosome ####
To simplify the dataset, isolate to a single chromosome and inspect across all populations.     

e.g., isolate to Chr18 (note: the source BCF file must be indexed for the following to work)       
```
# Subset the chromosome with bcftools 
bcftools view 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf --regions 18 -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.bcf 

```


To obtain only a region within the chromosome:      
```
# Index the VCF file
bcftools index 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz 

# Create a selected regions file
bcftools view 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz | grep -vE '^#' - | awk '$2 > 55000000 && $2 < 65000000 { print $1 "\t" $2 }' - > 02_input_data/selected_snp_region_file.txt

# Subset with selected regions file to produce compressed VCF file
bcftools view --regions-file 02_input_data/selected_snp_region_file.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz -Oz -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18_sel_region.vcf.gz

```

Plot a heatmap of the above selected file using the following RScript:    
`01_scripts/heatmap.R`








##### Extra info #####
To see raw genotypes of any specific locus, use the following code:    
```
# Replace <locus> with target position
bcftools view 02_input_data/*_miss0.15.vcf | grep -vE '^##' - | grep -E '^#|<locus>' - > 03_results/<locus>_genos.txt
```
Then open in a spreadsheet and transpose to get long form.    

