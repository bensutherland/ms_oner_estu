# Early Stuart sockeye salmon population genetics
Code to support analysis of Early Stuart sockeye metapopulation within the context of the upper Fraser River sockeye.     

### Whole-genome resequencing analysis ###
#### Requirements ####
- bcftools    
- amplitools 
- simple_pop_stats

Code is run from individual repositories, as indicated in the section.     
Clone `amplitools` and `simple_pop_stats`.     

#### Data sources ####
Download source VCF file [Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz](https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Christensen_et_al_2024/25705428).     
Put in `amplitools_*/02_input_data`.        

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

# Convert the BCF file to VCF file to be able to be read into R via vcfR    
bcftools view 14_extract_mhap/Oner.*_MAF0.05_5w50kb.bcf -Ov -o 14_extract_mhap/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_w_tags_MAF0.05_5w50kb.vcf

# Copy both of the VCF files back to simple_pop_stats
cp 14_extract_mhap/*.vcf ../simple_pop_stats/02_input_data/   
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


#### Population genetic analyses ####
Follow the instructions in the script `01_scripts/popgen_char_2025-02-27.R`.     


To see raw genotypes of any specific locus, use the following code:    
```
# Replace <locus> with target position
bcftools view 02_input_data/*_miss0.15.vcf | grep -vE '^##' - | grep -E '^#|<locus>' - > 03_results/<locus>_genos.txt
```
Then open in a spreadsheet and transpose to get long form.    

