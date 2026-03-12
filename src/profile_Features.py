import os
import sys
import argparse
from pathlib import Path

def main():

    parser = argparse.ArgumentParser(
        description="",
        usage="")

    parser.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS, Contigs, Coverage and Annotation nomenclature")
    
    parser.add_argument("--annotation",metavar="FILE",required=True,help="Annotation file containing your cds names and a feature (gene_name, COG, KO, EC_Number, etc. Must be a 2 column tab-separated file")
    parser.add_argument("--coverage",metavar="FILE",required=True,help="")
    parser.add_argument("--cov_mode",metavar="MODE",choices=['cds','contig'],default="cds",required=True,help="Format of your genomic coverage. It can be used as 'Coverage per Contigs' or 'Coverage per CDS' ")
    parser.add_argument("--feature_list",metavar="FILE",required=True,help="")
    parser.add_argument("--output_folder",metavar="DIR",required=True,help="")


    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()


    #coverage_file = open(sys.argv[1],"r")
    #tsc_name = str(sys.argv[3])
    #annotation = open(sys.argv[4],"r")
    #gene_list = open(sys.argv[5],"r")
    #tf_gen_res= str(sys.argv[6])
    #coverage_mode= str(sys.argv[7])

    #Coverages
    all_annot= dict()
    genes = dict()


    ##################################################################
    # Step 1. Read all genes of interest (must be on annotation)
    ##################################################################
    #print("Initiate gene dictionary")
    #The list should be in a format like:
    #gen1
    #gen2
    #gen3
    #gen4
    #The genes() dictionary should be:
    #gen1=0
    #gen2=0
    #gen3=0
    #The genes_array_tf should be:
    #[gen1,gen2,gen3]

    genes_array_tf = []
    with open(args.feature_list) as feature_list:
        for feature in feature_list:
            aux_g=feature.replace("\n","").lower()
            genes[aux_g]=0
            genes_array_tf.append(aux_g)


    genes_array_tf=list(genes_array_tf)
    genes_array_tf.sort()

    #print(genes_array_tf)

    ##################################################################
    # Step 2. Read all coverage of contigs.
    ##################################################################
    #print("Initiate coverage dictionary")
    #The list should be in a format like:
    #Contig_1[tab]Coverage_number
    #Contig_2[tab]Coverage_number
    #Contig_3[tab]Coverage_number
    #The result coverage_dict should be:
    #Contig_1=150
    #Contig_2=200
    #Contig_1=189

    coverage_dict=dict()
    with open(args.coverage) as coverage_file:
        for cover in coverage_file:
            aux=cover.split("\t")
            cov_value = aux[1].rstrip("\n")
            if("#" not in cover and args.cov_mode=="contig"):
                coverage_dict[aux[0]]=cov_value
            if("#" not in cover and args.cov_mode=="cds"):
                if(float(cov_value)>=0):
                    coverage_dict[aux[0]]=cov_value
                else:
                    coverage_dict[aux[0]]=0

    core=""
    depth=0
    gene_array=[]
    coverage_pairs=[]
    core_flag=False


    with open(args.annotation) as annotation:
        for lines in annotation:

            aux=lines.split("\t")
            cds_name=aux[0]
            responce_gene=aux[1].rstrip("\n").lower()
            scaffold_pair=cds_name.split("_")
            scaffold_pair="_".join(scaffold_pair[:-1])

            if("," in responce_gene):
                all_responce_genes=responce_gene.split(",")
                for any_responce in all_responce_genes:
                    any_responce=any_responce.rstrip("\n")
                    if(args.cov_mode=="contig"):
                        if(scaffold_pair in coverage_dict.keys() and any_responce in genes):
                            genes[any_responce]=genes[any_responce]+float(coverage_dict[scaffold_pair])
                        else:
                            continue
                    if(args.cov_mode=="cds"):
                        if(cds_name in coverage_dict.keys() and any_responce in genes):
                            genes[any_responce]=genes[any_responce]+float(coverage_dict[cds_name])
                            #print(cds_name, responce_gene, coverage_dict[cds_name])
                        else:
                            continue
            else:
                
                    if(args.cov_mode=="contig"):
                        if(scaffold_pair in coverage_dict.keys() and responce_gene in genes):
                            genes[responce_gene]=genes[responce_gene]+float(coverage_dict[scaffold_pair])
                        else:
                            continue
                    if(args.cov_mode=="cds"):
                        if(cds_name in coverage_dict.keys() and responce_gene in genes):
                            genes[responce_gene]=genes[responce_gene]+float(coverage_dict[cds_name])
                            print(cds_name, responce_gene, coverage_dict[cds_name])
                        else:
                            continue

    #################  S U M      O F       G E N E S          ################

    #print("#############################")
    #print("TOTAL SUM PER INSTANCES")


    all_data = genes.items()
    
    profiling_folder_only_feature = Path(args.output_folder+"/P_Features")
    profiling_folder_only_feature.mkdir(parents=True, exist_ok=True)
    output_folder_tf_gen = open(str(profiling_folder_only_feature)+"/"+args.sample_id+"_feature_profile.tsv","w")

    output_folder_tf_gen.write("Feature\t"+str(args.sample_id)+"\n")
    for item in all_data:
        output_folder_tf_gen.write(str(item[0])+"\t"+str(item[1])+"\n")

if __name__ == "__main__":
    main()


