# Early Stuart sockeye salmon population genetics
Code to support analysis of Early Stuart sockeye metapopulation

### 01. Whole-genome resequencing analysis ###
#### Requirements ####
- bcftools    
- simple_pop_stats
- amplitools 

Note: all code is run from the main directory of `simple_pop_stats`, unless otherwise specified.    
Clone both `simple_pop_stats` and `amplitools`.     

#### Data sources ####
Download source VCF file [Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz](https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Christensen_et_al_2024/25705428).     
Put in `02_input_data`.        

#### Reduce to only samples of interest ####
Keep only samples from populations of interest from the BCF file:    
```
bcftools query -l 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz | grep -E 'Bivouac|BowronR|Driftwood|Dust|Felix|Kuzkwa|MiddleR|Nadina|Paula|Pinchi|Stellako|Tachie|Takla_S' - > 02_input_data/samples_to_retain.txt

bcftools view -S 02_input_data/samples_to_retain.txt 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr.vcf.gz -Ob -o 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf    

# Copy file to amplitools
cp 02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained.bcf ../amplitools/14_extract_mhap/
```

#### Filter target BCF file ####
This section is executed from the amplitools directory.     

```
# Filter the BCF file as per amplitools standards
01_scripts/filter_bcf.sh

# Follow amplitools README to add all annotations and to filter based on MAF (MAF 0.05).    
# Caution: intermediate files will fill up storage rapidly! (update once resolved in amplitools)

# Generate an LD filtered dataset for some analyses:    
bcftools +prune -w 50kb -m 0.5 -Ob -o 14_extract_mhap/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_w_tags_MAF0.05_5w50kb.bcf 14_extract_mhap/Oner*_w_tags_MAF0.05.bcf`     

# Convert the BCF file to VCF file to be able to be read into R via vcfR    
bcftools view 14_extract_mhap/Oner.*_MAF0.05_5w50kb.bcf -Ov -o 14_extract_mhap/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_w_tags_MAF0.05_5w50kb.vcf

bcftools view 14_extract_mhap/Oner.*_MAF0.05.bcf -Ov -o 14_extract_mhap/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_w_tags_MAF0.05.vcf

# Copy both of the VCF files back to simple_pop_stats
cp 14_extract_mhap/*.vcf ../simple_pop_stats/02_input_data/   
```

#### Population genetic analyses ####
Follow the instructions in the script `01_scripts/popgen_char_2025-02-27.R`.     


To see raw genotypes of any specific locus, use the following code:    
```
# Replace <locus> with target position
bcftools view 02_input_data/*_miss0.15.vcf | grep -vE '^##' - | grep -E '^#|<locus>' - > 03_results/<locus>_genos.txt
```
Then open in a spreadsheet and transpose to get long form.    

