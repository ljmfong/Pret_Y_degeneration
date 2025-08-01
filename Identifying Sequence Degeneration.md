# Identifying Sequence Degeneration:

Here are the commands used to identify the sequence degeneration of _P. reticulata_ following the attainment of Dthe NA sequences from the respective wild populations (Almeida et al. 2021) and lab population (Lin et al. 2023). 
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

We used **snpEFF** to predict the effect of a SNP as "low", "moderate", or "high" impact on gene function (
