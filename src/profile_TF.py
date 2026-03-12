import os
import sys
import numpy as np
import argparse
from pathlib import Path

def main ():

    parser = argparse.ArgumentParser(
        description="",
        usage="")

    parser.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS, Contigs, Coverage and Annotation nomenclature")
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


    # Step 1. Read all coverage of contigs.
    #print("Initiate coverage dictionary")

    coverage_dict=dict()
    with open(args.coverage) as coverage_file:
        for cover in coverage_file:
            aux=cover.split("\t")
            if("#" not in cover and args.cov_mode=="contig"):
                coverage_dict[aux[0]]=aux[1].rstrip("\n")
            if("#" not in cover and args.cov_mode=="cds"):
                coverage_dict[aux[0]]=aux[1].rstrip("\n")

    # Step 2. Read the TFMs list
    #print("Read TFMs dictionary nomenclature")

    TFMs_dict=dict()
    TFMs_array= []
    
    with open(args.tf_list) as tfm_list:
        for tfs in tfm_list:
            aux=tfs.split("\t")
            #print(tfs)
            TFMs_dict[aux[0]]=aux[1].rstrip("\n") # Nomenclature
            #print(aux[0])
            TFMs_array.append(aux[1].rstrip("\n")) # To create the array op TFMs positions


    TFMs_array = list(set(list(TFMs_array)))
    TFMs_array=sorted(TFMs_array,key=str.lower)

    #print (TFMs_array)
    #TFMs_array.sort()

    counting_TFMs=dict()

    for t_array in TFMs_array:
        counting_TFMs[t_array]=0
    #       print(ble)
    #print(TFMs_array)
    #print(TFMs_dict


    # Step 3. Read the TFMs per Operon file, with all the gene info.
    #print("Generating the TFMs sums")


    all_operon_data=[]
    operon_dict_cds = dict()
    operon_dict_tf = dict()


    with open(args.prr_results) as prr_tfs:
        for prr_line in prr_tfs:
            aux_array=""
            aux=prr_line.split(",")
            tf_name=aux[0].replace(".txt","")
            if(tf_name in TFMs_dict):
                simple_TFMs=TFMs_dict[tf_name]
                aux_opr=aux[1]
                aux_contig=aux[2].split("_")        
                aux_contig="_".join(aux_contig[:-1])
                evalue=float(aux[-1].rstrip("\n"))
                #print(evalue)
                if(args.cov_mode=="contig"):
                    if (evalue<=float(args.cutoff)):  
                        aux_array=simple_TFMs+","+aux_opr+","+aux_contig
                        #print(aux_array)
                        all_operon_data.append(aux_array)
                if(args.cov_mode=="cds"): 
                    if (evalue>=float(args.cutoff)):  
                        if(aux_opr not in operon_dict_cds.keys()):
                            operon_dict_cds.setdefault(aux_opr,[])
                            operon_dict_cds[aux_opr].append(aux[2])
                            
                            operon_dict_tf.setdefault(aux_opr,[])
                            operon_dict_tf[aux_opr].append(simple_TFMs)
                        else:
                            operon_dict_cds[aux_opr].append(aux[2])
                            operon_dict_tf[aux_opr].append(simple_TFMs)
                                         

    #print(len(all_operon_data))
    all_operon_data = list(set(list(all_operon_data)))
    #print(len(all_operon_data))


    if(args.cov_mode=="contig"):
        for triplet in all_operon_data:
            aux=triplet.split(",")
            #print(triplet)
            #print(aux[0])
            contig=aux[2]
            counting_TFMs[aux[0]]=float(counting_TFMs[aux[0]])+float(coverage_dict[contig])
            #print(counting_TFMs[aux[0]])
            #print(str(coverage_dict[contig]))
            #geneTF[aux[2]]=aux[0].replace(".txt","")       
            #print(len(aux))

    if(args.cov_mode=="cds"):
        for key, values in operon_dict_tf.items():
            tf_list=list(set(values))
            uniq_genes = list(set(operon_dict_cds[key]))
            coverage_operon=[]
            for any_gene in uniq_genes:
                coverage_operon.append(float(coverage_dict[any_gene]))
            coverage_operon_mean = np.mean(coverage_operon) 
            for any_tf in tf_list:
                counting_TFMs[any_tf]=float(counting_TFMs[any_tf]) + float(coverage_operon_mean)
            #print(key)
            #print(uniq_genes)
            #print(tf_list)
            #print(coverage_operon)
            #print(coverage_operon_mean)


    profiling_folder_only_tf = Path(args.output+"/P_TF")
    profiling_folder_only_tf.mkdir(parents=True, exist_ok=True)
    output_tf_gen = open(str(profiling_folder_only_tf)+"/"+args.sample_id+
 ".tsv","w")
    output_tf_gen.write("TF_Group\t"+str(args.sample_id)+"\n")
    for keys in counting_TFMs.keys():
        output_tf_gen.write(keys+"\t"+str(counting_TFMs[keys])+"\n")

if __name__ == "__main__":
    main()



