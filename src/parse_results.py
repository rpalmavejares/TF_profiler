import sys
import argparse

def main ():

    parser = argparse.ArgumentParser(
        description="",
        usage="")

    parser.add_argument("--prr_fasta",metavar="",help="")    
    parser.add_argument("--motif_alignment",metavar="",help="")
    parser.add_argument("--sample_id",metavar="",help="")
    parser.add_argument("--motif_db",metavar="",help="")
    parser.add_argument("--output",metavar="",help="")

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()    

    operon_cds_dict=dict()
    operon_contig_dict=dict()
    
    with open(args.prr_fasta) as prr_fasta_file:
        for operons in prr_fasta_file:
            if(">" in operons):
                aux=operons.split(" ")
                operon_name=aux[0].replace(">","")
                contig_name=aux[1].replace("Replicon=","")
                genes=aux[5].replace("Genes=","")
                operon_contig_dict[operon_name]=contig_name
                operon_cds_dict[operon_name]=genes.replace("\n","")

    TF_flag=False
    Best=False
    Register=False

    all_data=[]

    motif_db_dir=(args.motif_db).rstrip("/")
   
    outfile = open(args.output,"a")
 
    with open(args.motif_alignment) as motif_alignment_file:
        for mast in motif_alignment_file:
            if("# mast -norc -nostatus -hit_list -best" in mast):
                if(len(all_data)>0):
                    TFBS_taxa=mast.split(" ")
                    TFBS_taxa=TFBS_taxa[6].replace(motif_db_dir+"/","")
                    for data in all_data:
                        outfile.write(TFBS_taxa+","+data+"\n")	
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
                if("," in genes):
                    each_gene=genes.split(",")
                    for each in each_gene:
                        all_data.append(aux[0]+","+each+","+",".join(mast_stats[3:10]))
                else:
                    all_data.append(aux[0]+","+genes+","+",".join(mast_stats[3:10]))

        if(len(all_data)>0):
            for data in all_data:
                outfile.write(TFBS_taxa+","+data,"\n")

if __name__ == "__main__":
    main()










	
