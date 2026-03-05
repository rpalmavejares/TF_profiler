import os
import sys
import numpy as np

#TSC258 TARA_B110000211_G
#TSC000  TARA_X000000368_G

coverage_file = open(sys.argv[1],"r")
station = str(sys.argv[2]) 
tsc_name = str(sys.argv[3])
annotation = open(sys.argv[4],"r")
gene_list = open(sys.argv[5],"r")
#logs_file= open("TFMs_constanza/outlayers/"+sys.argv[3]+".txt","w")
tfm_list= open(sys.argv[6],"r")
operon_tfs= open(sys.argv[7],"r")
tf_tf_res= str(sys.argv[8])
tf_gen_res= str(sys.argv[9])
e_cutoff= int(sys.argv[10])
coverage_mode = str(sys.argv[11])


tf_tf_results = tf_tf_res+tsc_name+".mtx"
tf_gen_results = tf_gen_res+tsc_name+".mtx"


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
for gen in gene_list:
    aux_g=gen.replace("\n","").lower()
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
for cover in coverage_file:
    aux=cover.split("\t")
    if("#" not in cover and coverage_mode=="contig"):
        coverage_dict[aux[0]]=aux[6].replace("\n","")
    if("#" not in cover and coverage_mode=="cds"):
        if(float(aux[5])>=50):
            coverage_dict[aux[0]]=aux[6].replace("\n","")
        else:
            coverage_dict[aux[0]]=0



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
for tfs in tfm_list:
    aux=tfs.split("\t")
    TFMs_dict[aux[0]]=aux[1].replace("\n","") # Nomenclature
    TFMs_array.append(aux[1].replace("\n","")) # To create the array op TFMs positions


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
for operon_line in operon_tfs:
    aux=operon_line.split(",")
    geneTF.setdefault(aux[2],[])
    operonTF.setdefault(aux[1],[])
    gene_operonTF.setdefault(aux[2])
    gene_operonTF[aux[2]]=aux[1]
    operon_dict_cds.setdefault(aux[1],[])
    operon_dict_cds[aux[1]].append(float(coverage_dict[aux[2]]))
    #print(aux[0])
    aux_TFMs=aux[0].replace(".txt","")
    evalue=aux[-1].split("e-")
    #print (evalue[1])

    if (int(evalue[1].replace("\n",""))>=e_cutoff): 
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
    #print(str(item[0])+"\t"+str(item[1]))
    uniqss = list(set(list(item[1])))
    operonTFMs[item[0]]=uniqss
    operon_once[item[0]]=uniqss
    #print(item[0],uniqs) 


#print("######################")
#print(geneTFMs) # PAIRS OF TARA GENE AND TF
#print("######################")
#print(operonTFMs) # PAIRS OF OPERON AND TF
#print("######################")
#print(gene_operonTF) # PAIRS OF OPERON AND GENE
#print("######################")


#print (geneTF["TARA_B110000211_G_scaffold46761_1_gene50144"])






##################################################################
# Step 6. Read all clusters from CD-HIT
#print("Initiate finder of Cluster"+"\t"+str(tsc_name))
##################################################################


cluster_list=[]
all_cluster=[]
core=""
depth=0
tf_per_gene = [0] * len(TFMs_array)
gene_array=[]


#print(gene_array)





##################################################################
# Step 7. Create matrix for the TFMs and Gene Vector
##################################################################
gene_array = [[0 for x in range(len(TFMs_array))] for y in range(len(genes_array_tf)) ]
TF_array = [[0 for x in range(len(TFMs_array))] for y in range(len(genes_array_tf)) ]

#gene_array[0][694]=100
coverage_pairs=[]
core_flag=False

#print(genes)
#print(len(genes))



for lines in annotation:

    aux=lines.split("\t")
    cds_name=aux[0]
    responce_gene=aux[1].rstrip("\n").lower()
    scaffold_pair=cds_name.split("_")
    scaffold_pair="_".join(scaffold_pair[:-1])

    if(coverage_mode=="contig"):
        if(cds_name in geneTFMs and responce_gene in genes):
            genes[responce_gene]=genes[responce_gene]+float(coverage_dict[scaffold_pair])
            
            #if(responce_gene=="agua"):
            #print(subcluster_gene+"\t"+str(float(coverage_dict[scaffold_pair[0]])))
            #count=0
            #print(str(coverage_dict[scaffold_pair[0]]))
            operon_value=gene_operonTF[cds_name]
            for TFMs in geneTFMs[cds_name]:
                #if(responce_gene=="agua" and TFMs=="FUR"):
                #   print(subcluster_gene)
                operon_value=gene_operonTF[cds_name]
                TF_ID=TFMs_array.index(TFMs)
                GENE_ID=genes_array_tf.index(responce_gene)
                gene_array[GENE_ID][TF_ID]= gene_array[GENE_ID][TF_ID] + float(coverage_dict[scaffold_pair])
            if(operon_value in operon_once.keys()):
                for TFMs in geneTFMs[cds_name]:
                    #if(responce_gene=="agua" and TFMs=="FUR"):
                    #   print(subcluster_gene)
                    operon_value=gene_operonTF[cds_name]    
                    TF_ID=TFMs_array.index(TFMs)
                    GENE_ID=genes_array_tf.index(responce_gene)
                    TF_array[GENE_ID][TF_ID]=TF_array[GENE_ID][TF_ID] + float(coverage_dict[scaffold_pair])
                operon_once.pop(operon_value)
        else:
            continue
    if(coverage_mode=="cds"):

        if(cds_name in geneTFMs and responce_gene in genes):
            
            operon_value=gene_operonTF[cds_name]
            new_avg_coverage=np.mean(np.array(operon_dict_cds[operon_value]))
            #if(responce_gene=="agua"):
            #print(subcluster_gene+"\t"+str(float(coverage_dict[scaffold_pair[0]])))
            #count=0
            #print(str(coverage_dict[scaffold_pair[0]]))
            
            genes[responce_gene]=genes[responce_gene]+float(new_avg_coverage)
            operon_value=gene_operonTF[cds_name]
            for TFMs in geneTFMs[cds_name]:
                #if(responce_gene=="agua" and TFMs=="FUR"):
                #   print(subcluster_gene)
                operon_value=gene_operonTF[cds_name]
                TF_ID=TFMs_array.index(TFMs)
                GENE_ID=genes_array_tf.index(responce_gene)
                gene_array[GENE_ID][TF_ID]= gene_array[GENE_ID][TF_ID] + float(new_avg_coverage)
            if(operon_value in operon_once.keys()):
                for TFMs in geneTFMs[cds_name]:
                    #if(responce_gene=="agua" and TFMs=="FUR"):
                    #   print(subcluster_gene)
                    operon_value=gene_operonTF[cds_name]
                    TF_ID=TFMs_array.index(TFMs)
                    GENE_ID=genes_array_tf.index(responce_gene)
                    TF_array[GENE_ID][TF_ID]=TF_array[GENE_ID][TF_ID] + float(new_avg_coverage)
                operon_once.pop(operon_value)
        else:
            continue




#print(genes_array_tf)
final_count_TFMs=[0] * len(TFMs_array)


############     V E C T O R     T F M S      BY     G E N    #################

#print("#############################")
#print("TFMs VECTOR PER INSTANCE")

print("Instance\t"+str("\t".join(TFMs_array)))
ocurrences=len(gene_array[0])
for data in range(len(gene_array)):
    if(gene_array[data].count(0)!=ocurrences):
        print(str(genes_array_tf[data])+"\t"+"\t".join(str(v) for v in gene_array[data]))      # VECTOR DE TFMS POR GEN
    for TFMs_pos in range(len(gene_array[data])):   
        final_count_TFMs[TFMs_pos]=float(final_count_TFMs[TFMs_pos]+TF_array[data][TFMs_pos])


#print("#############################")
#print("SUM OF TFMs")


###########          S U M       O F       T F M S        ############

output_tf_tf = open(tf_tf_results,"w")
output_tf_tf.write("Instance\t"+str(tsc_name)+"\n")
for values in range(len(final_count_TFMs)):
    output_tf_tf.write(str(TFMs_array[values])+"\t"+str(final_count_TFMs[values])+"\n")

output_tf_tf.close()
#print(final_count_TFMs)

#all_data = genes.items()


#################  S U M      O F       G E N E S          ################

#print("#############################")
#print("TOTAL SUM PER INSTANCES")


all_data = genes.items()
output_tf_gen = open(tf_gen_results,"w")

output_tf_gen.write("Instance\t"+str(tsc_name)+"\n")
for item in all_data:    # RESULTADO DE SUMA DE COVERTURAS POR GEN
    output_tf_gen.write(str(item[0])+"\t"+str(item[1])+"\n")




