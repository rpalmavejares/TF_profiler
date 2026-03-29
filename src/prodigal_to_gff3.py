import sys
import argparse

def main():


    parser = argparse.ArgumentParser(
        description="",
        usage="")

    parser.add_argument("--input",metavar="FILE", required=True,help="Prodigal cds file containing all sequences. Nucleotide or Amoniacidic formats can be used.")
    parser.add_argument("--output",metavar="FILE", required=True,help="Output file for the GFF file")

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    outfile= open(args.output,"w")
    outfile.write("##gff-version 3\n")
    outfile.write("# Companion tool for TF_profiler\n")
    outfile.write("# Simplified notation for CDS sequences\n")

    with open(args.input,"r") as prodigal_file:
        for lines in prodigal_file:
            if(">" in lines):
                aux= lines.rstrip("\n").split(" # ")
                parent_contig = aux[0].lstrip(">")
                parent_contig = "_".join((parent_contig.split("_"))[:-1] )
                cds_name = aux[0].lstrip(">")
                start = aux[1]
                stop = aux[2]
                if(aux[3]=="1"):
                    strand = "+"
                if(aux[3]=="-1"):
                    strand = "-"
                outfile.write(str(parent_contig)+"\tTF_profiler\tCDS\t"+str(start)+"\t"+str(stop)+"\t.\t"+str(strand)+"\t.\tID="+str(cds_name)+"\n")




if __name__ == "__main__":
    main()
