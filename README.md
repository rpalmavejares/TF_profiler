# TF profiler
TF profiler is a tool created to search TF motifs into genomic material (MetaG, MAGs, etc) given a simplified regulation model. The tool will calculate de abundace of Motifs that regulate CDS inside marked operons on the genomic source.

## TF Profiler Pipeline and What it does.

The pipeline is based on 3 Stages, which are described on the image below.

As a short explanation:

### TF Profile: Stage 1.

We obtain all the CDS from contigs in Metagenomes or Single Genomes (MAGs) provided by the user.
These contigs must not contain "N" or "unknown" base pairs, so is imperative to not use scaffolds with uncertain nucleotides.

For each all CDS marked on the contigs we build a operons model based on their Intergenic distance. If 2 CDS are separated by less than 50 bp, we consider them part of the same operon.
This calculation is made by "contigs strand" which means that the model makes a 2 pass caltulation, 3' to 5' and 5' to 3'.

Once all operons are generated, we mark the PRR or Potential Regulatory Region, defined as the nucleotide segment in the intergenic region in front of the operon. This process is tied to each strand, same as before.
The size of the PRR is shown in the model as 300bp upstream of the first CDS in the operon, and 30bp downstrean of the first CDS in the operon.
This PRR is also generated with operon with only 1 CDS of size.

### TF Profile: Stage 2.

### TF Profile: Stage 3.

<br>

![A descriptive alt text for your image](imgs/TF_profiler_Full_Model.jpg)


### TF Profile: Regulatory Matrix.

### Running the Pipeline.

