import os
import sys
import numpy as np
import argparse
from pathlib import Path
import numpy as np

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
    counting_Features = dict()

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
    #The feature_array_tf should be:
    #[gen1,gen2,gen3]

    feature_array_tf = []
    with open(args.feature_list) as gene_list:
        for gen in gene_list:
            aux_g=gen.replace("\n","").lower()
            counting_Features[aux_g]=0
            feature_array_tf.append(aux_g)


    feature_array_tf=list(feature_array_tf)
    feature_array_tf.sort()


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

    cds_feature= dict()

    with open(args.annotation) as annotation:
        for annot in annotation:
            aux_annot=annot.split("\t")
            cds_id=aux_annot[0]
            feature_name=aux_annot[1].rstrip("\n")
            if(args.sample_id in cds_id and feature_name not in ("_", "", " ", "-", "*", ".")):
                cds_feature[cds_id]=feature_name.lower()
            
    #print(cds_feature)

    operon_and_cds = dict()
    operon_and_tf = dict()
    operon_and_test = dict()
    count_call=0
    count_pre_call=0
    with open(args.prr_results) as prr_tfs:
        for prr_line in prr_tfs:
            aux_array=""
            aux=prr_line.split(",")
            tf_name=aux[0].replace(".txt","")
            if(tf_name in TFMs_dict):
                simple_TFMs=TFMs_dict[tf_name]
                aux_opr=aux[1]
                aux_cds=aux[2]
                evalue=float(aux[-1].rstrip("\n"))
                aux_eva = aux[-1].rstrip("\n")
                aux_eva = aux_eva.split("e-")
                if(int(aux_eva[1])>=6):
                    count_pre_call+=1
                if(evalue<=float(args.cutoff)):
                    count_call+=1
                    if(aux_opr not in operon_and_cds): # works for both dicts
                        operon_and_cds.setdefault(aux_opr,[])
                        operon_and_tf.setdefault(aux_opr,[])
                        operon_and_cds[aux_opr].append(aux_cds)
                        operon_and_tf[aux_opr].append(simple_TFMs)
                    else:
                        operon_and_cds[aux_opr].append(aux_cds)
                        operon_and_tf[aux_opr].append(simple_TFMs)
    
        
    counting_TFMs=dict()
    for t_array in TFMs_array:
        counting_TFMs[t_array]=0
    counter = 0

    TF_Feature_vectors = np.zeros((len(counting_Features),len(counting_TFMs))) 

    for key, values in operon_and_tf.items():
        all_tfs_list = list(set(values))
        all_cds_list = list(set(operon_and_cds[key]))
        cds_check=False
        if(args.cov_mode == "contig"):
            parent_contig=all_cds_list[0].split("_")
            parent_contig="_".join(parent_contig[:-1])
            for any_cds in all_cds_list:
                if(any_cds in cds_feature.keys()):
                    if cds_feature[any_cds] in counting_Features.keys():
                        if(parent_contig in coverage_dict.keys()):
                            #print("Check-True",key, parent_contig, all_cds_list)
                            cds_check=True
                        else:
                            print(f'An error has ocurred {parent_contig} not in coverage files')
            if(cds_check and len(all_cds_list)>=1):
                # Countig all TFs in an Operon
                for any_tfs in all_tfs_list:
                    counting_TFMs[any_tfs]= float(counting_TFMs[any_tfs]) + float(coverage_dict[parent_contig])
                cds_check=False
                
                # Countig all cds Features in an Operon
                
                for any_cds in all_cds_list:
                    if any_cds in cds_feature.keys():
                        feature_name = cds_feature[any_cds]
                        counting_Features[feature_name]= float(counting_Features[feature_name]) + float(coverage_dict[parent_contig])
                        for any_ffs_m in all_tfs_list:
                            TFMs_index = TFMs_array.index(any_ffs_m)
                            Feature_index = feature_array_tf.index(feature_name)
                            TF_Feature_vectors[Feature_index][TFMs_index]=float(TF_Feature_vectors[Feature_index][TFMs_index]) + float(coverage_dict[parent_contig])
                

                #print(key, all_tfs_list)
                #print(key, all_cds_list)
                #break
        
        if(args.cov_mode == "cds"):
            operon_cds_coverages=[]
            for any_cds in all_cds_list:
                if(any_cds in cds_features.keys()):
                    if(cds_feature[any_cds] in counting_Features.keys()):
                        if(any_cds in coverage_dict.keys()):
                            operon_cds_coverages.append(float(coverage_dict[any_cds]))
                            cds_check=True
                        else:
                            print(f'An error has ocurred {any_cds} not in coverage files')
            if(cds_check and len(all_cds_list)>=1):
                mean_coverage = np.mean(operon_cds_coverages) 
                # Countig all TFs in an Operon
                for any_tfs in all_tfs_list:
                    counting_TFMs[any_tfs]= float(counting_TFMs[any_tfs]) + float(mean_coverage)
            
                # Countig all cds Features in an Operon
                
                for any_cds in all_cds_list:
                    if any_cds in cds_feature.keys():
                        feature_name = cds_feature[any_cds]
                        counting_Features[feature_name]= float(counting_Feature[feature_name]) + float(coverage_dict[any_cds])
                        for any_ffs_m in all_tfs_list:
                            TFMs_index = TFMs_array.index(any_ffs_m)
                            Feature_index = feature_array_tf.index(feature_name)
                            TF_Feature_vectors[Feature_index][TFMs_index]=float(TF_Feature_vectors[Feature_index][TFMs_index]) + float(coverage_dict[any_cds])




    #for key, values in counting_TFMs.items():
    #    print(f'{key}\t{values}')

    #for key, values in counting_Features.items():
    #    print(f'{key}\t{values}')


    print("Feature/TF\t"+str("\t".join(TFMs_array)))
    for rows in range(len(TF_Feature_vectors)):
        if(np.sum(TF_Feature_vectors[rows])>0):
            #print(np.sum(rows))
            Feature_index = feature_array_tf[rows]
            rows_parsed = "\t".join(map(str,TF_Feature_vectors[rows]))
            print(f'{Feature_index}\t{rows_parsed}')




        #print(key, all_tfs_list)
        #print(key, all_cds_list)
        #for tfs_s in all_tfs_single:
        #    print(key)        

   
    print("count_pre_calls",count_pre_call)
    print("count_calls",count_call)
 
    sys.exit(0)

 

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
            output_tf_feature_vectors.write(str(feature_array_tf[data])+"\t"+"\t".join(str(v) for v in gene_array[data]))
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


    all_data = counting_Features.items()
    output_tf_gen = open(tf_gen_results,"w")

    profiling_folder_tf_feature_feat = Path(args.output+"/P_TF-Features/By_Feature")
    profiling_folder_tf_feature_feat.mkdir(parents=True, exist_ok=True)
    output_tf_feature_feat = open(str(profiling_folder_tf_feature_feat)+"/"+args.sample_id+"_tf_profile.tsv","w")
    
    output_tf_feature_tfms.write("TF_Group\t"+str(args.sample_id)+"\n")
    for item in all_data:    # RESULTADO DE SUMA DE COVERTURAS POR GEN
        output_tf_feature_feat.write(str(item[0])+"\t"+str(item[1])+"\n")

if __name__ == "__main__":
    main()


