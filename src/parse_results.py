import sys

fasta_operon_file = open(sys.argv[1],"r")
mast_alignment_file = open(sys.argv[2],"r")


operon_cds_dict=dict()
operon_contig_dict=dict()




for operons in fasta_operon_file:

	if(">" in operons):
		aux=operons.split(" ")
		operon_name=aux[0].replace(">","")
		contig_name=aux[1].replace("Replicon=","")
		genes=aux[5].replace("Genes=","")
		operon_contig_dict[operon_name]=contig_name
		operon_cds_dict[operon_name]=genes.replace("\n","")
		#print(operon_name,contig_name,genes)

TF_flag=False
Best=False
Register=False

all_data=[]

for mast in mast_alignment_file:
	if("# mast -norc -nostatus -hit_list -best Regprecise_TF_DB/" in mast):
		if(len(all_data)>0):
			TFBS_taxa=mast.split(" ")
			TFBS_taxa=TFBS_taxa[6].replace("Regprecise_TF_DB/","")
			#print(TFBS_taxa,all_data)
			for data in all_data:
				print(TFBS_taxa+","+data)	
			all_data=[] 
		continue
	if("# Best single (non-overlapping)" in mast):
		Best=True
		continue
	if("# sequence_name" in mast):
		Register=True
		continue
	if(Best and Register):
		aux=mast.split(" ")
		genes=operon_cds_dict[aux[0]]
		mast_stats=[]
		for values in aux:
			if(values!=""):
				mast_stats.append(values.replace("\n",""))
		#print(mast_stats)	
		if("," in genes):
			each_gene=genes.split(",")
			for each in each_gene:
				#print(aux[0]+","+each+","+",".join(mast_stats[1:10]))
				all_data.append(aux[0]+","+each+","+",".join(mast_stats[3:10]))
		else:
			#print(aux[0]+","+genes+","+",".join(mast_stats[1:10]))
			all_data.append(aux[0]+","+genes+","+",".join(mast_stats[3:10]))

if(len(all_data)>0):
	for data in all_data:
		print(TFBS_taxa+","+data)











	
