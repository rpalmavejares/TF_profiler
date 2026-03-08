import os
import sys
import numpy as np
import argparse
from pathlib import Path

#coverage_file = open(sys.argv[1],"r")
#station = str(sys.argv[2]) 
#annotation = open(sys.argv[4],"r")
#gene_list = open(sys.argv[5],"r")
#logs_file= open("TFMs_constanza/outlayers/"+sys.argv[3]+".txt","w")
#tfm_list= open(sys.argv[6],"r")
#operon_tfs= open(sys.argv[7],"r")
#tf_tf_res= str(sys.argv[8])
#tf_gen_res= str(sys.argv[9])
#e_cutoff= int(sys.argv[10])
#coverage_mode = str(sys.argv[11])

def main ():


    parser = argparse.ArgumentParser(
        description="",
        usage="")

    parser.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS, Contigs, Coverage and Annotation nomenclature")
    parser.add_argument("--annotation",metavar="FILE",required=True,help="Annotation file containing your cds names and a feature (gene_name, COG, KO, EC_Number, etc. Must be a 2 column tab-separated file")
    parser.add_argument("--feature_list",metavar="FILE",required=True,help="")
    parser.add_argument("--coverage",metavar="FILE",required=True,help="")
    parser.add_argument("--cov_mode",metavar="MODE",choices=['cds','contig'],default="cds",required=True,help="Format of your genomic coverage. It can be used as 'Coverage per Contigs' or 'Coverage per CDS' ")
    parser.add_argument("--cutoff",metavar="FLOAT",type=float,default=1e6,required=True,help="E-value cutoff for mapped motifs agains PRR")
    parser.add_argument("--tf_list",metavar="FILE",required=True,help="")
    parser.add_argument("--prr_results",metavar="FILE",required=True,help="")
    parser.add_argument("--output",metavar="DIR",required=True,help="")

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

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
    with open(args.feature_list) as gene_list:
        for gen in gene_list:
            aux_g=gen.replace("\n","").lower()
            genes[aux_g]=0
            genes_array_tf.append(aux_g)


    genes_array_tf=list(genes_array_tf)
    genes_array_tf.sort()


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
            if("#" not in cover and args.cov_mode=="contig"):
                coverage_dict[aux[0]]=aux[1].rstrip("\n")
            if("#" not in cover and args.cov_mode=="cds"):
                coverage_dict[aux[0]]=aux[1].rstrip("\n")


    ##################################################################
    # Step 3. Read the TFMs list
    ##################################################################
    #print("Read TFMs dictionary nomenclature")
    #The list should be in a format like:
    #TFM1_long_name[tab]TFM1_short_name 
    #TFM2_long_name[tab]TFM2_short_name
    #TFM3_long_name[tab]TFM3_short_name
    #The result dictionary should be:
    #TFM1_long_name=TFM1_short_name
    #TFM2_long_name=TFM2_short_name
    #TFM3_long_name=TFM3_short_name
    #The TFMs_array should be like:
    #sorted[TFM1,TFM2,TFM3,...]


    TFMs_dict=dict()
    TFMs_array= []
    with open(args.tf_list) as tfm_list:
        for tfs in tfm_list:
            aux=tfs.split("\t")
            TFMs_dict[aux[0]]=aux[1].rstrip("\n") # Nomenclature
            TFMs_array.append(aux[1].rstrip("\n")) # To create the array op TFMs positions


    TFMs_array = list(set(list(TFMs_array)))
    TFMs_array.sort()


    ##################################################################
    # Step 4. Read the TFMs per Operon file, with all the gene info.
    ##################################################################
    # print("Reading the TFMs pero operon file, from mast results")
    #The list should be in a format like:
    #TFM_subclass.txt,operon_id,tara_gene_name,strand,start,stop,score,evalue
    #AcnR_Corynebacteriaceae.txt,O_50_3,TARA_B110001469_G_scaffold124594_2_gene138560,+1,11,26,-178.00,4.23e-05
    #The geneTF dictionary should be: 
    #tara_gene_name[TFM1,TFM5]
    #tara_gene_name[TFM7,TFM1,TFM2]
    #The operonTF dictionary should be:
    #operon_id[TFM1,TFM5]
    #operon_id[TFM7,TFM1,TFM2]
    #The gene_operonTF dictionary should be:
    #tara_genen_name[operon_id]
    #tara_genen_name[operon_id]


    geneTF=dict()
    operonTF=dict()
    gene_operonTF=dict()
    operon_dict_cds=dict()
    operon_dict_tf = dict()   
 
    with open(args.prr_results) as prr_tfs:
        for prr_line in prr_tfs:
            aux=prr_line.split(",")
            geneTF.setdefault(aux[2],[])
            operonTF.setdefault(aux[1],[])
            gene_operonTF.setdefault(aux[2])
            gene_operonTF[aux[2]]=aux[1]
            operon_dict_cds.setdefault(aux[1],[])
            print(aux[0],aux[1],aux[2])
            operon_dict_cds[aux[1]].append(float(coverage_dict[aux[2]]))
            aux_TFMs=aux[0].replace(".txt","")
            evalue=float(aux[-1].rstrip("\n"))


            if (evalue>=float(args.cutoff)): 
                if(aux_TFMs in TFMs_dict):
                    simple_TFMs=TFMs_dict[aux_TFMs]
                    geneTF[aux[2]].append(simple_TFMs)
                    operonTF[aux[1]].append(simple_TFMs)
            else:
                continue
            #geneTF[aux[2]]=aux[0].replace(".txt","")   
            #print(len(aux))


    ##################################################################
    # Step 5. Read all annotation data from eggnogmapper
    ##################################################################
    #print("Initiate annotation dictionary")
    #The list should be in a format like:
    #query[tab]seed_ortholog[tab]evalue[tab]score[tab]eggnog_ogs[tab]max_annot_lvl[tab]cog[tab]description[tab]preferred_name[tab]gos[tab]ec
    #or#
    #query[tab]prefered_name
    #The all_annot dictionary should be:
    #query(tara_gene)=preferred_name
    #query(tara_gene)=preferred_name
    #query(tara_gene)=preferred_name
    #The annotation_excluded dictionary should be:
    #If the evalue leaves a preferred_name out, then is stored as:
    #query(tara_gene)=preferred_name
    #query(tara_gene)=preferred_name
    #print(annotation_excluded)

    all_data = geneTF.items()
    geneTFMs=dict()

    #print("geneTF : ",geneTF)
    for item in all_data:
            #print(str(item[0])+"\t"+str(item[1]))
        uniqs = list(set(list(item[1])))
        geneTFMs[item[0]]=uniqs
        if(len(uniqs)==0):
            geneTFMs.pop(item[0])
        #print(item[0],uniqs)

    all_data_operon = operonTF.items()
    operonTFMs=dict()
    operon_once=dict()
    for item in all_data_operon:
        uniqss = list(set(list(item[1])))
        operonTFMs[item[0]]=uniqss
        operon_once[item[0]]=uniqss

    #print("######################")
    #print(geneTFMs) # PAIRS OF TARA GENE AND TF
    #print("######################")
    #print(operonTFMs) # PAIRS OF OPERON AND TF
    #print("######################")
    #print(gene_operonTF) # PAIRS OF OPERON AND GENE
    #print("######################")

    core=""
    depth=0
    tf_per_gene = [0] * len(TFMs_array)
    gene_array=[]



    ##################################################################
    # Step 7. Create matrix for the TFMs and Gene Vector
    ##################################################################
    gene_array = [[0 for x in range(len(TFMs_array))] for y in range(len(genes_array_tf)) ]
    TF_array = [[0 for x in range(len(TFMs_array))] for y in range(len(genes_array_tf)) ]

    coverage_pairs=[]
    core_flag=False

    with open(args.annotation) as annotation:
        for lines in annotation:

            aux=lines.split("\t")
            cds_name=aux[0]
            responce_gene=aux[1].rstrip("\n").lower()
            scaffold_pair=cds_name.split("_")
            scaffold_pair="_".join(scaffold_pair[:-1])

            if(args.cov_mode=="contig"):
                if(cds_name in geneTFMs and responce_gene in genes):
                    genes[responce_gene]=genes[responce_gene]+float(coverage_dict[scaffold_pair])
                    
                    operon_value=gene_operonTF[cds_name]
                    for TFMs in geneTFMs[cds_name]:
                        operon_value=gene_operonTF[cds_name]
                        TF_ID=TFMs_array.index(TFMs)
                        GENE_ID=genes_array_tf.index(responce_gene)
                        gene_array[GENE_ID][TF_ID]= gene_array[GENE_ID][TF_ID] + float(coverage_dict[scaffold_pair])
                    if(operon_value in operon_once.keys()):
                        for TFMs in geneTFMs[cds_name]:
                            operon_value=gene_operonTF[cds_name]    
                            TF_ID=TFMs_array.index(TFMs)
                            GENE_ID=genes_array_tf.index(responce_gene)
                            TF_array[GENE_ID][TF_ID]=TF_array[GENE_ID][TF_ID] + float(coverage_dict[scaffold_pair])
                        operon_once.pop(operon_value)
                else:
                    continue
            if(args.cov_mode=="cds"):

                if(cds_name in geneTFMs and responce_gene in genes):
                    
                    operon_value=gene_operonTF[cds_name]
                    new_avg_coverage=np.mean(np.array(operon_dict_cds[operon_value]))
                    genes[responce_gene]=genes[responce_gene]+float(new_avg_coverage)
                    operon_value=gene_operonTF[cds_name]
                    for TFMs in geneTFMs[cds_name]:
                        operon_value=gene_operonTF[cds_name]
                        TF_ID=TFMs_array.index(TFMs)
                        GENE_ID=genes_array_tf.index(responce_gene)
                        gene_array[GENE_ID][TF_ID]= gene_array[GENE_ID][TF_ID] + float(new_avg_coverage)
                    if(operon_value in operon_once.keys()):
                        for TFMs in geneTFMs[cds_name]:
                            operon_value=gene_operonTF[cds_name]
                            TF_ID=TFMs_array.index(TFMs)
                            GENE_ID=genes_array_tf.index(responce_gene)
                            TF_array[GENE_ID][TF_ID]=TF_array[GENE_ID][TF_ID] + float(new_avg_coverage)
                        operon_once.pop(operon_value)
                else:
                    continue

    final_count_TFMs=[0] * len(TFMs_array)


    ############     V E C T O R     T F M S      BY     G E N    ############

    #print("#############################")
    #print("TFMs VECTOR PER FEATURE")


    profiling_folder_tf_feature_vectors = Path(args.output+"/P_TF-Features/Vectors")
    profiling_folder_tf_feature_vectors.mkdir(parents=True, exist_ok=True)
    output_tf_feature_vectors = open(str(profiling_folder_tf_feature_vectors)+"/"+args.sample_id+"_tf-feature_vectors_profile.tsv","w")
    output_tf_feature_vectors.write("TF_Group\t"+str("\t".join(TFMs_array)))
    
    ocurrences=len(gene_array[0])
    for data in range(len(gene_array)):
        if(gene_array[data].count(0)!=ocurrences):
            output_tf_feature_vectors.write(str(genes_array_tf[data])+"\t"+"\t".join(str(v) for v in gene_array[data]))
            output_tf_feature_vectors.write("\n")      # VECTOR DE TFMS POR GEN
        for TFMs_pos in range(len(gene_array[data])):   
            final_count_TFMs[TFMs_pos]=float(final_count_TFMs[TFMs_pos]+TF_array[data][TFMs_pos])

    output_tf-feature_vectors.close()

    #print("#############################")
    #print("SUM OF TFMs")


    ###########          S U M       O F       T F M S        ############

    profiling_folder_tf_feature_tfms = Path(args.output+"/P_TF-Features/By_TF")
    profiling_folder_tf_feature_tfms.mkdir(parents=True, exist_ok=True)
    output_tf_feature_tfms = open(str(profiling_folder_tf_feature_tfms)+"/"+args.sample_id+"_tf_profile.tsv","w")
    output_tf_feature_tfms.write("TF_Group\t"+str(args.sample_id)+"\n")

    for values in range(len(final_count_TFMs)):
        output_tf_feature_tfms.write(str(TFMs_array[values])+"\t"+str(final_count_TFMs[values])+"\n")
    output_tf_feature_tfms.close()



    #################  S U M      O F       G E N E S          ################

    #print("#############################")
    #print("TOTAL SUM PER INSTANCES")


    all_data = genes.items()
    output_tf_gen = open(tf_gen_results,"w")

    profiling_folder_tf_feature_feat = Path(args.output+"/P_TF-Features/By_Feature")
    profiling_folder_tf_feature_feat.mkdir(parents=True, exist_ok=True)
    output_tf_feature_feat = open(str(profiling_folder_tf_feature_feat)+"/"+args.sample_id+"_tf_profile.tsv","w")
    
    output_tf_feature_tfms.write("TF_Group\t"+str(args.sample_id)+"\n")
    for item in all_data:    # RESULTADO DE SUMA DE COVERTURAS POR GEN
        output_tf_feature_feat.write(str(item[0])+"\t"+str(item[1])+"\n")

if __name__ == "__main__":
    main()


