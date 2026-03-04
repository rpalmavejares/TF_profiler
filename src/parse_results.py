import sys
import argparse
from pathlib import Path

def main ():

    parser = argparse.ArgumentParser(
        description="Parser the results of motifs mapping and convert them to a simplified csv format",
        usage="python %(prog)s --prr_fasta <> ")

    parser.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser.add_argument("--prr_fasta",metavar="FILE",type=Path,help="Path to file containing PRR sequences")
    parser.add_argument("--motif_alignment",metavar="FILE",help="Path to file containing the MAST aligment output")
    parser.add_argument("--output",metavar="FILE",help="Output file of parsed results")

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

   
    outfile = open(args.output,"a")
 
    with open(args.motif_alignment) as motif_alignment_file:
        for mast in motif_alignment_file:
            if("# mast -norc -nostatus -hit_list -best" in mast):
                if(len(all_data)>0):
                    TFBS_taxa=mast.split(" ")
                    TFBS_taxa = TFBS_taxa[6].split("/")
                    TFBS_taxa = TFBS_taxa[-1]
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










	
