import os
import sys
import numpy as np

#TSC258 TARA_B110000211_G

coverage_file = open(sys.argv[1],"r")
tfm_list = open(sys.argv[2],"r") 
operon_tfs = open(sys.argv[3],"r")
station= str(sys.argv[4])
eva_cutoff= str(sys.argv[5])
coverage_mode = str(sys.argv[6])

# Step 1. Read all coverage of contigs.
#print("Initiate coverage dictionary")

coverage_dict=dict()
for cover in coverage_file:
    aux=cover.split("\t")
    if("#" not in cover and coverage_mode=="contig"):
        coverage_dict[aux[0]]=aux[6].replace("\n","")
    if("#" not in cover and coverage_mode=="cds"):
        if(float(aux[5])>=50):
            coverage_dict[aux[0]]=aux[6].replace("\n","")
        else:
            coverage_dict[aux[0]]=0
# Step 2. Read the TFMs list
#print("Read TFMs dictionary nomenclature")

TFMs_dict=dict()
TFMs_array= []
for tfs in tfm_list:
    aux=tfs.split("\t")
    #print(tfs)
    TFMs_dict[aux[0]]=aux[1].replace("\n","") # Nomenclature
    #print(aux[0])
    TFMs_array.append(aux[1].replace("\n","")) # To create the array op TFMs positions


TFMs_array = list(set(list(TFMs_array)))
TFMs_array=sorted(TFMs_array,key=str.lower)

#print (TFMs_array)
#TFMs_array.sort()

counting_TFMs=dict()

for ble in TFMs_array:
    counting_TFMs[ble]=0
#       print(ble)
#print(TFMs_array)
#print(TFMs_dict


# Step 3. Read the TFMs per Operon file, with all the gene info.
#print("Generating the TFMs sums")


all_operon_data=[]
operon_dict_cds = dict()
operon_dict_tf = dict()


for operon_line in operon_tfs:
    aux_array=""
    aux=operon_line.split(",")
    tf_name=aux[0].replace(".txt","")
    if(tf_name in TFMs_dict):
        simple_TFMs=TFMs_dict[tf_name]
        aux_opr=aux[1]
        aux_contig=aux[2].split("_")        
        aux_contig="_".join(aux_contig[:-1])
        aux_evalue=aux[-1]
        evalue=aux[-1].split("e-")
        #print (evalue[1])
        if(coverage_mode=="contig"):
            if (int(evalue[1].replace("\n",""))>=int(eva_cutoff)):  
                aux_array=simple_TFMs+","+aux_opr+","+aux_contig
                #print(aux_array)
                all_operon_data.append(aux_array)
        if(coverage_mode=="cds"): 
            if (int(evalue[1].replace("\n",""))>=int(eva_cutoff)):  
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


if(coverage_mode=="contig"):
    for triplet in all_operon_data:
        aux=triplet.split(",")
        #print(triplet)
        #print(aux[0])
        contig=aux[2]
        counting_TFMs[aux[0]]=float(counting_TFMs[aux[0]])+float(coverage_dict[contig])
            #geneTF[aux[2]]=aux[0].replace(".txt","")       
            #print(len(aux))

if(coverage_mode=="cds"):
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



print ("TFM\t"+station) 
for keys in counting_TFMs.keys():
    print (keys+"\t"+str(counting_TFMs[keys]))








