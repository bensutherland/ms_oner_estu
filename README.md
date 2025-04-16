# Early Stuart sockeye salmon population genetics
Code to support analysis of Early Stuart sockeye metapopulation within the context of the upper Fraser River sockeye.     

### Whole-genome resequencing analysis ###
#### Requirements ####
- bcftools    
- amplitools 
- simple_pop_stats

Commands are run from individual repositories, as indicated in the section.     
Clone `amplitools` and `simple_pop_stats`.     

#### Data sources ####
Download source VCF file [Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz](https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Christensen_et_al_2024/25705428) and put it in `amplitools_*/02_input_data`.        

#### Subset samples to target region ####
Keep only samples from populations of interest from the BCF file:    
```
# Start by choosing scope, either Option A or Option B, as listed below

# Option A: Stuart and Nadina/Stellako only
bcftools query -l 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz | grep -E 'Bivouac|BowronR|Driftwood|Dust|Felix|Kuzkwa|MiddleR|Nadina|Paula|Pinchi|Stellako|Tachie|Takla_S' - > 02_input_data/samples_to_retain.txt

# Option B: All upper Fraser River
bcftools query -l 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz | grep -E 'Bivouac|Driftwood|Dust|Felix|Paula|Takla|Kuzkwa|Middle|Pinchi|Tachie|Stellako|Nadina|Bowron|Horsefly|BlueLead|McKinley|Mitchell|horsefly|Wasko|Quesnel|Chilko|Taseko' - | grep -vE '\_K\_' - > 02_input_data/samples_to_retain.txt

# Use sample file to subset the VCF file to only these samples 
bcftools view -S 02_input_data/samples_to_retain.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf    

```

#### Filter variants in subset BCF file ####
This section is executed from the amplitools directory.     

Filter using `01_scripts/filter_bcf.sh` by updating the source folder, BCF file name, number cores, other parameters as needed.    

Additional filters:       
```
# Create filtered dataset for population genetic uses
# MAF filter, first add tags
bcftools +fill-tags 02_input_data/Oner.*_maxDP10000_minGQ20_miss0.15.bcf -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q20_avgDP7_biallele_minDP7_maxDP10000_minGQ20_miss0.15_w_tags.bcf

# Then filter on tags for MAF
bcftools view -i 'MAF > 0.05' 02_input_data/*_w_tags.bcf -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q20_avgDP7_biallele_minDP7_maxDP10000_minGQ20_miss0.15_w_tags_MAF0.05.bcf

# Generate an LD filtered, compressed VCF file for population genetic analyses     
bcftools +prune -w 50kb -m 0.5 -Oz -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q20_avgDP7_biallele_minDP7_maxDP10000_minGQ20_miss0.15_w_tags_MAF0.05_5w50kb.vcf.gz 02_input_data/*_MAF0.05.bcf

# Copy the VCF file to simple_pop_stats
cp 02_input_data/*5w50kb.vcf.gz ../simple_pop_stats/02_input_data/   

# Change directory to simple_pop_stats
cd ../simple_pop_stats_v.0.2/
```

Additional dataset: create another VCF file specific to only the four collections from EStu that have the most samples, then conduct MAF and LD filters. This will be used for population genetic analysis within EStu. Do this from within `simple_pop_stats`.        
```
# Select only the populations of interest
bcftools query -l 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15.vcf | grep -E 'Driftwood|Dust|Bivouac|Paula' - > 02_input_data/samples_to_retain_estu_limited.txt

bcftools view -S 02_input_data/samples_to_retain_estu_limited.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15.vcf -Ov -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited.vcf

# Fill tags and filter on MAF
bcftools +fill-tags 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited.vcf -Ov -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags.vcf

bcftools view -i 'MAF > 0.01' 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags.vcf -Ov -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags_MAF0.01.vcf

# Filter on LD
bcftools +prune -w 50kb -m 0.5 -Ov -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags_MAF0.01_5w50kb.vcf 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags_MAF0.01.vcf

```
...this will be used as an input in the below RScript.    


Additional dataset: Chilcotin, selection dataset
```
# Create sample list from the all loci file
bcftools query -l 02_input_data/*_w_tags.bcf | grep -E 'Chilko|Taseko' - > 02_input_data/samples_to_retain_chilcotin.txt

# Subset
bcftools view -S 02_input_data/samples_to_retain_chilcotin.txt 02_input_data/*_miss0.15_w_tags.bcf -Ob -o 02_input_data/chilcotin.bcf

# Fill tags again to recalculate and filter on MAF
bcftools +fill-tags 02_input_data/chilcotin.bcf -Ob -o 02_input_data/chilcotin_w_tags.bcf
bcftools view -i 'MAF > 0.01' 02_input_data/chilcotin_w_tags.bcf -Oz -o 02_input_data/chilcotin_w_tags_maf0.01.vcf.gz

# Clean space
rm 02_input_data/chilcotin.bcf 02_input_data/chilcotin_w_tags.bcf


# Copy to simple_pop_stats

# This will be used by 01_scripts/chilcotin.R
```

Additional dataset: Quesnel/Horsefly, selection dataset
```
# Create sample list from the all loci file
bcftools query -l 02_input_data/*_w_tags.bcf | grep -E 'Horsefly|BlueLead|McKinley|Mitchell|horsefly|Wasko|Quesnel' - > 02_input_data/samples_to_retain_quesnel.txt

# Subset
bcftools view -S 02_input_data/samples_to_retain_quesnel.txt 02_input_data/*_miss0.15_w_tags.bcf -Ob -o 02_input_data/quesnel.bcf

# Fill tags again to recalculate and filter on MAF
bcftools +fill-tags 02_input_data/quesnel.bcf -Ob -o 02_input_data/quesnel_w_tags.bcf
bcftools view -i 'MAF > 0.01' 02_input_data/quesnel_w_tags.bcf -Oz -o 02_input_data/quesnel_w_tags_maf0.01.vcf.gz
bcftools view -i 'MAF > 0.05' 02_input_data/quesnel_w_tags.bcf -Oz -o 02_input_data/quesnel_w_tags_maf0.05.vcf.gz

# For testing, produce a LD-filtered
bcftools +prune -w 50kb -m 0.5 -Oz -o 02_input_data/quesnel_w_tags_maf0.05_5w50kb.vcf.gz 02_input_data/quesnel_w_tags_maf0.05.vcf.gz

# Clean space
rm 02_input_data/quesnel.bcf 02_input_data/quesnel_w_tags.bcf

# Move outputs to simple_pop_stats

# This will be used by 01_scripts/quesnel.R
```

Additional dataset: Stuart, selection dataset      
```
# Create sample list from the all loci file
bcftools query -l 02_input_data/*_w_tags.bcf | grep -E 'Bivouac|Driftwood|Dust|Felix|Paula|Takla|Kuzkwa|Middle|Pinchi|Tachie' - > 02_input_data/samples_to_retain_stuart.txt

# Subset
bcftools view -S 02_input_data/samples_to_retain_stuart.txt 02_input_data/*_miss0.15_w_tags.bcf -Ob -o 02_input_data/stuart.bcf

# Selected chromosome only
bcftools index 02_input_data/stuart.bcf
bcftools view 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf --regions 18 -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.bcf 

# MAF filter (light)
bcftools +fill-tags 02_input_data/stuart_chr18.bcf -Ob -o 02_input_data/stuart_chr18_w_tags.bcf
bcftools view -i 'MAF > 0.01' 02_input_data/stuart_chr18_w_tags.bcf -Oz -o 02_input_data/stuart_chr18_w_tags_maf0.01.vcf.gz

# Move outputs to simple_pop_stats
```






#### All populations, all loci, one chromosome ####
To simplify the dataset, isolate to a single chromosome and inspect across all populations.     

e.g., isolate to Chr18 (note: the source BCF file must be indexed for the following to work)       
```
# Subset the chromosome with bcftools 
bcftools view 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf --regions 18 -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.bcf 

```


To obtain only a region within the chromosome:      
```
# Create a selected regions file
bcftools view 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz | grep -vE '^#' - | awk '$2 > 55000000 && $2 < 65000000 { print $1 "\t" $2 }' - > 02_input_data/selected_snp_region_file.txt

# Subset with selected regions file to produce compressed VCF file
bcftools view --regions-file 02_input_data/selected_snp_region_file.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz -Oz -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18_sel_region.vcf.gz

```







#### Population genetic analyses ####
Follow the instructions in the script `01_scripts/popgen_char_2025-02-27.R`.     


To see raw genotypes of any specific locus, use the following code:    
```
# Replace <locus> with target position
bcftools view 02_input_data/*_miss0.15.vcf | grep -vE '^##' - | grep -E '^#|<locus>' - > 03_results/<locus>_genos.txt
```
Then open in a spreadsheet and transpose to get long form.    

