# Identifying Sequence Degeneration:

Here are the commands used to identify the sequence degeneration of _P. reticulata_ following the attainment of Dthe NA sequences from the respective wild populations (Almeida et al. 2021) and lab population (Lin et al. 2023). 
----------------------------------------------------------

## Alignment to female reference genome

We used the Ensembl female _P. reticulata_ reference genome and aligned the paired-end DNA sequences using **BWA** and sorted using **SAMTools**. Here is an example command:

    bwa mem -t 12 Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa sample_forward_R1.fastq.gz sample_reverse_R2.fastq.gz > sample_id.sam
    samtools sort -o sample_id.sorted.sam sample_id.sam

After all the samples were mapped to the reference genome and sorted
