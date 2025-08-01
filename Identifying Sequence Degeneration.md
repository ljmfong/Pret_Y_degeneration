# Identifying Sequence Degeneration:

Here are the commands used to identify the sequence degeneration of _P. reticulata_ following the attainment of Dthe NA sequences from the respective wild populations (Almeida et al. 2021) and lab population (Lin et al. 2023). 
----------------------------------------------------------

## Alignment to female reference genome

We used the Ensembl female _P. reticulata_ reference genome and aligned the DNA sequences using *BWA*. Here is an example command:

  /Linux/bin/bwa mem -t 12 ~/assemblies/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa /MankFlex/ljmfong/Yuying_ret_family/NS.1357.002.IDT_i7_9---IDT_i5_9.P433_R1.fastq.gz /MankFlex/ljmfong/Yuying_ret_family/NS.1357.002.IDT_i7_9---IDT_i5_9.P433_R2.fastq.gz > /MankFlex/ljmfong/yuying_pret/Fam1_Mother_P433.sam
  /Linux/samtools-1.9/bin/samtools sort -o /MankFlex/ljmfong/yuying_pret/sorted_sam_females/Fam1_Mother_P433.sorted.sam /MankFlex/ljmfong/yuying_pret/Fam1_Mother_P433.sam
