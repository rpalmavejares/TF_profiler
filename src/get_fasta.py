import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import argparse



def main():

    parser = argparse.ArgumentParser(
        description="Extract the Potential Regulatory Region from the Simplified Operon Model output",
        usage="python %(prog)s --sample_id <ID> --operon_model <FILE> --output_folder <DIR>")

    parser.add_argument("--sample_id",metavar="ID",help="Sample ID or Name for creating the temporary files")
    parser.add_argument("--operon_model",metavar="FILE",help="Path to file containing the prr model")
    parser.add_argument("--output_folder",metavar="DIR",help="Output folder for all intermediate files")    

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    output_folder=args.output_folder.rstrip("/")
    output_fasta = open(output_folder+"/"+args.sample_id+"_prr.fasta","w")

    with open(args.output_folder+"/"+args.operon_model) as prr_file:
        all_records = []
        for lines in prr_file:
            genes=list()
            aux = lines.replace("\n","")
            aux=aux.split(",")
            replicon = aux[0]
            operon_name = aux[1]
            strand = aux[2]
            sequence = aux[3]
            start = aux[4]
            stop = aux[5]
            genes = aux[6:]
            if(len(sequence)!=0):
                description="Replicon="+str(replicon)+" Strand="+str(strand)+" Start="+str(start)+" Stop="+str(stop)+" Genes="+str(",".join(genes))

                #print(operon_name)
                record = SeqRecord(Seq(sequence), id=operon_name, description=description)
                all_records.append(record)

        SeqIO.write(all_records,output_fasta,"fasta")
    output_fasta.close()
if __name__ == "__main__":
    main()
