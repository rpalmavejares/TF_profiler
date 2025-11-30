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

Finally, with the use of MAST (MEME suite) we map a set of Transcriptipn Factor (TF) motifs coming from the database RegPrecise. https://regprecise.lbl.gov/collections_tf.jsp to all our PRR created on the previous steps.
These motifs are categorised in 88 macro groups or TF.
Once the mapping of all motifs is finished, a final file wit core information is generated called "motif profile".

### TF Profile: Stage 2 and Stage 3

These Stages are run togheter, but they are easier to explain on 2 separated ways.

#### Stage 2:

This step takes the file created in the previous step (motif profile) and parses it together with 2 new input. 
* CDS Features File:
  <br>
  A tab separated file containing the CDS IDs and a genomic description. This could be gene_names, COG categories, KO pathways, EC numbers, OG, etc.

```
CDS_IDs                              GENE_NAME
sample_1_contig_Nº_cds_Nº                    -
sample_1_contig_Nº_cds_Nº+1                liga
sample_1_contig_Nº_cds_Nº+2                 pkg
sample_1_contig_Nº_cds_Nº+3                 fur
sample_1_contig_Nº_cds_Nº+4                  -

```



  * Contigs / CDS Coverage


<br>

![A descriptive alt text for your image](imgs/TF_profiler_Full_Model.jpg)



### Running the Pipeline.

