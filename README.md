# TF profiler
TF profiler is a tool created to search TF motifs into genomic material (MetaG, MAGs, etc) given a simplified regulation model. The tool will calculate de abundace of Motifs that regulate CDS inside marked operons on the genomic source.

## TF Profiler Pipeline.

The pipeline is based on 3 Stages, which are described on the image below.

<br>

![A descriptive alt text for your image](imgs/TF_profiler_Full_Model.jpg)

## Pipeline description:

### TF Profile: Stage 1.

The first step is to obtain all CDS from contigs in Metagenomes or Single Genomes (MAGs) provided by the user.
These contigs must not contain "N" or "unknown" base pairs, so is imperative to not use scaffolds with uncertain nucleotides.

For each all CDS marked on the contigs we build a operons model based on their Intergenic distance. If 2 CDS are separated by less than 50 bp, we consider them part of the same operon.
This calculation is made by "contigs strand" which means that the model makes a 2 pass caltulation, 3' to 5' and 5' to 3'.

Once all operons are generated, we mark the PRR or Potential Regulatory Region, defined as the nucleotide segment in the intergenic region in front of the operon. This process is tied to each strand, same as before.
The size of the PRR is shown in the model as 300bp upstream of the first CDS in the operon, and 30bp downstrean of the first CDS in the operon.
This PRR is also generated with operon with only 1 CDS of size.

Finally, with the use of MAST (MEME suite) we map a set of Transcriptipn Factor (TF) motifs coming from the database RegPrecise. https://regprecise.lbl.gov/collections_tf.jsp to all our PRR created on the previous steps.
These motifs are categorised in 88 macro groups or TF.
Once the mapping of all motifs is finished, a final file wit core information is generated called "motif profile".

### TF Profile: Stage 2.

The second step in the pipeline is to filter for (or select) the TF motifs that bind to all PRRs based on a given E-value cutoff from the MAST mappings in the previous step. The highest E-value reported by MAST is $10^{-4}$; however, those values represent less than 0.01% of the total hits, as most fall within the $10^{-5}$ cutoff. On the other hand, even though the lowest value can theoretically be 0, our observations show a minimum value of $10^{-20}$. We strongly recommend using an E-value cutoff between $10^{-5}$ and $10^{-7}$, as lower values are too strict for filtering purposes. Once an E-value cutoff is chosen, the TF motif hits and their respective PRRs are merged with the CDS Annotation File, which contains descriptive features. This data integration follows this structure:

1) A **Transcription Factor motifs** mapped to a **PRR**
2) A **PRR** that belongs to an **OPERON**
3) An **OPERON** that contains a/many **CDS**
4) A **CDS** that contains a descriptive **FEATURE**, such as gene_name, KO pathways, EC Numbers, COG, etc.

Finally, after all the data is combined, a **TARGET File** is used to select relevant features. This TARGET File must contain features already present in the CDS Annotation File.

### TF Profile: Stage 3.

The last step of the pipeline consists of adding the Coverage File (for Contigs or CDS) to calculate the total abundance of each TF motif and each TARGET feature in our sample.

*  For Metagenomes, we recommend the usage of Contig coverages.
*  For Metatranscriptomes, we recommend the usage of CDS coverages.

Regarding the information structure, we add a final field: 

5) A **CDS** that contains a **COVERAGE** value (quantified by contig or CDS).

Within this structure, we have implemented the following rules:

*  If a **TF motif** appears two or more times in any given **PRR**, it is only **counted once**.
*  If an **OPERON** contains two or more CDS, the coverage of the **TFs** in the **PRR** is also **counted once**.

Those 2 rules try to minimize the impact of:

*  **Database Bias:** More popular TFs often have more motif entries.
*  **Operon Size Bias:** Small and large operons should be treated as a single regulated unit, rather than being counted once per CDS




##  TF Profile Inputs.

For running all 3 Stages we require the following files:

*  A Fasta file containing the genomic data.
*  A PTT file containing the positions of all CDS.
*  An Annotation file containing the features of each CDS.
*  A Coverage file containing a metric of abundance by Contig or CDS.
*  A TARGET file containing the features to keep in our data structure.


At this stage we require 2 different files
<br>

* Fasta File with Contigs
  
```
>Sample_001_contig_0001
GTACGTGACTCGTGACTTACGTCAGGGCGGCGGGGATTACTTACTTA............
>Sample_001_contig_0002
AACTTTAGCTTGCGGATCTATTCGTACCACATGCTATCTGATTTCTT............
>Sample_001_contig_0003
TCTTTCGAQTATTATCGGATCTTACGGTGTATGACATACACACTTAG............
```

* Protein Table File (PTT file) with CDS positions by contig.

```
Metagenome assembly, contig: Sample_001_contig_0001
2 proteins
Replicon  Location  Strand  Length  PID  Gene  Synonym Code  COG  Product
Sample_001_contig_0001  77..472  -  395  Sample_001_contig_0001_1  -  Sample_001_contig_0001_1  -  -  hypothetical protein
Sample_001_contig_0001  474..986  -  512  Sample_001_contig_0001_2  -  Sample_001_contig_0001_2  -  -  hypothetical protein
```

#### Create a PTT File from a FAA Protein File (Prodigal).

If you already have a .faa file with the CDS from Prodigal, you can create a PTT file with the provided script.

Example Prodigal .faa file:

```
>Sample_001_contig_0001_1  # 77 # 472 # -1 # 
MYELKEYLNAINVSKESLLDSEDEMWEKKYAPFIVNKCVAPFPDTILLVNEVNQYHHLDK
KLQFDFLLNSLRTRKRYTPWLKAKKLKNLEYVKEYYGYNNEKAKAALDILDDEQISAIKI
KLNKGGRNGRN*
>Sample_001_contig_0001_2  # 474 # 986 # -1 # 
MIDIYDNVLEPHLAELIDLNLKQQTWKYDYHSQQGTPNKHWHVFCGHNPWEVTSKDYEWL
MPIWDTALAKYNFKEKYNVSEFKRLYLNAHTHGIEPHMHIDDGDFTMMYYPRLDWKMDWG
GGTVVDGQLVQNIGNRLIVFPAYAPHQAQPVSRQCYDLRTVVVFKTWVDK*
```

With this file you can run the ptt_from_faa.py from the script folder:

```
ptt_from_faa.py Prodigal_File.faa > PTT_File.ptt
```

### TF Profile: Stage 2 and Stage 3

These Stages are run togheter, but they are easier to explain on 2 separated ways.

#### Stage 2:

This step takes the file created in the previous step (motif profile) and parses it together with 2 new input. 
* CDS Features File:
  <br>
  A tab separated file containing the CDS IDs and a genomic description. This could be gene_names, COG categories, KO pathways, EC numbers, OG, etc.

```
CDS_IDs                              GENE_NAME
sample-1_contig-Nº_cds-Nº                    -
sample-1_contig-Nº_cds-Nº+1                liga
sample-1_contig-Nº_cds-Nº+2                 pkg
sample-1_contig-Nº_cds-Nº+3                 fur
sample-1_contig-Nº_cds-Nº+4                  -
```




  * Contigs / CDS Coverage





### Running the Pipeline.

