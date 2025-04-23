
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

