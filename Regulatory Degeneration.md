# Identifying Regulatory Degeneration:

Here are the commands used to identify the regulatory degeneration of _P. reticulata_ following the attainment of scRNA-seq data from Darolti & Mank (2023). 
----------------------------------------------------------

## Alignment to female reference genome

We used the Ensembl female _P. reticulata_ reference genome and aligned the paired-end DNA sequences using **BWA** and sorted using **SAMTools**. Here is an example command:

    bwa mem -t 12 Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa sample_forward_R1.fastq.gz sample_reverse_R2.fastq.gz > sample_id.sam
    samtools sort -o sample_id.sorted.sam sample_id.sam

After all the samples were mapped to the reference genome and sorted, I made a text file containing a list of the sorted bam files and their full paths. Then I used **BCFTools** to create my VCF and indexed the VCF.

        bcftools mpileup --threads 12 -Ou -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -f Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa -b Pret_bamlist.txt | bcftools call -mv -Oz -f GQ,GP -o pret_individ.vcf.gz
        tabix pret_individ.vcf.gz	
        
        
For the wild populations, the individuals were batched into 4 regions then merged together:

        ref=ref_f.fa
        region_name=$1 #batch calling, I seperate it into 4 regions.
        bcftools=~/bin/bcftools-1.16/bcftools
        ${bcftools} mpileup --threads 20 -Ou -q 20 -Q 20 -R region${region_name}.txt -b bamlist.txt -a FORMAT/AD,FORMAT/DP -f $ref  | ${bcftools} call -mv -Oz -a GQ,GP -o all.fem.LGs${region_name}.vcf.gz
        tabix all.fem.all_LGs.indels.vcf.gz

Then all files were merged together:

        bcftools merge -m all --threads 12 all.fem.all_LGs.indels.vcf.gz pret_individ.vcf.gz -o all_unfiltered_pret_snps.vcf -O z
        bcftools view all_unfiltered_pret_snps.vcf | bgzip -c  > all_unfiltered_pret_snps.vcf.gz
        tabix all_unfiltered_pret_snps.vcf.gz		

## Filtering and calling effect of SNPs and indels

We used **snpEFF** to predict the effect of a SNP as "low", "moderate", or "high" impact on gene function (details can be found **here**), and filtreing was done using **VCFTools**. The VCFs for SNPs and INDELS were done separately to be assessed separately. 

        /Linux/jdk-17.0.5/bin/java -jar ~/snpEff/snpEff.jar download -v Guppy_female_1.0_MT.99
        /Linux/jdk-17.0.5/bin/java -Xmx8g -jar ~/snpEff/snpEff.jar -v -stats snpEff_unfiltered_genomic_all_SNPs.html Guppy_female_1.0_MT.99 all_unfiltered_pret_snps.vcf.gz > genomic_SNPs_unfiltered.vcf
        vcftools --vcf genomic_SNPs_unfiltered.vcf --recode --recode-INFO-all --remove-indels --min-meanDP 7 --max-meanDP 112 --min-alleles 2 --max-alleles 3 --maf 0.05 --minQ 30 --max-missing 0.9 --out genomic_SNPs_filtered_no_indels
        vcftools --vcf genomic_SNPs_unfiltered.vcf --recode --recode-INFO-all --keep-only-indels --min-meanDP 7 --max-meanDP 112  --min-alleles 2 --max-alleles 3 --maf 0.05 --minQ 30 --max-missing 0.9 --out genomic_SNPs_filtered_indels

I separated out the different impacts of the SNPs based on the annotated effect by snpEFF using **snpSift**, examples provided below:
Low impact:

        /Linux/jdk-17.0.5/bin/java -Xmx8g -jar ~/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'coding_sequence_variant'" ~/new_Ret_SNPs/snpEff_results/genomic_SNPs_filtered_no_indels.recode.vcf > coding_sequence_variant_variants.vcf	

Moderate impact:

        /Linux/jdk-17.0.5/bin/java -Xmx8g -jar ~/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'missense_variant'" ~/new_Ret_SNPs/snpEff_results/genomic_SNPs_filtered_no_indels.recode.vcf > missense_variants.vcf

High impact:

        /Linux/jdk-17.0.5/bin/java -Xmx8g -jar ~/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'stop_gained'" ~/new_Ret_SNPs/snpEff_results/genomic_SNPs_filtered_no_indels.recode.vcf > stop_gained_variants.vcf


Files were then concatenated and prepared to run using the Rscript:

        cat high_impact_VCFs > high_impact_variants.vcf
        sed -i '/^##/d' high_impact_variants.vcf
        wc -l high_impact_variants.vcf

## Checking for duplicate genes:

I used VCFTools to identify duplicates following the pipeline by **Lin et al. (2022)**:

        vcftools --gzvcf all_unfiltered_pret_snps.vcf.gz --site-depth --keep male_individs.txt --out male_unfilt_markdupe_filtered	
        vcftools --gzvcf all_unfiltered_pret_snps.vcf.gz --site-depth --keep female_individs.txt --out female_unfilt_markdupe_filtered	
        python gene_coverage.py gene_boundary.bed male_unfilt_markdupe_filtered.ldepth female_unfilt_markdupe_filtered.ldepth unfilt_markdupe_MFDepth.txt			
        



