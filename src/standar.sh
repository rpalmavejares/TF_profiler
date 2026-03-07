#!/bin/bash
#SBATCH --job-name=Chile
#SBATCH --partition=main
#SBATCH --output=log_chile/out.%a.%N.%j.out
#SBATCH --error=log_chile/out.%a.%N.%j.err
##SBATCH --mail-user=rpalmavejares@gmail.com
##SBATCH --mail-type=ALL
#SBATCH -c 1 # number of cores
#SBATCH --mem=50G

##module load java/jdk1.8u102

#ml MEME/5.0.4-Perl-5.28.1-Python-2.7.15

intergenic_start=300
intergenic_stop=30
default_intergenic_zone=$2
gene_distance=$2


#output_folder="intergenic_150-10_TFMs" 
#output_folder="intergenic_150-10_TF" 
#output_folder="intergenic_300-30_TFMs" 
#output_folder="intergenic_contigs" 
#output_folder="intergenic_400-50_TF" 
ptt_folder="Tara_Chile_all_ptt"
fasta_folder="Tara_Chile_all_fasta"
output_folder="Testing_Pipeline" 

mkdir $output_folder


#python new_intergenic_pipeline/get_intergenics.py $ptt_folder/$1.ptt $fasta_folder/$1_500bp.fasta $intergenic_start $intergenic_stop $default_intergenic_zone $gene_distance > $output_folder/$1.intergenics_$gene_distance.tsv

#python new_intergenic_pipeline/get_fasta.py $output_folder/$1.intergenics_$gene_distance.tsv $output_folder/$1.operons_$gene_distance.fasta

#perl new_intergenic_pipeline/run_mast.pl mast Regprecise_TF_DB $output_folder/$1_mast_search_$gene_distance Regprecise_TF_DB/motifs.list $output_folder/$1.operons_$gene_distance.fasta

#python new_intergenic_pipeline/parse_results.py $output_folder/$1.operons_$gene_distance.fasta $output_folder/$1_mast_search_$gene_distance.txt > $output_folder/$1.operon_$gene_distance.npal 

#echo "succesull $1"

#python tf_profiler.py calculate --sample_id S9_Z05 --cds_pos ../Tara_Chile_all_ptt/S9_Z05.ptt --assembly ../Tara_Chile_all_fasta/S9_Z05_500bp.fasta --output_folder ../map_out --motif_list src/motifs.list

#python profile_Features.py --sample_id S9_Z05 --annotation ../../map_out/annotation_simple.tsv --coverage ../../map_out/coverage_simple.tsv --cov_mode contig --targets ../../map_out/all_gene_list.tsv --output M4
