import sys
import argparse



def main():


    parser = argparse.ArgumentParser(
        description="",
        usage="")
    
    parser.add_argument("--input",metavar="File",help="Temporal file with merged profiles")
    parser.add_argument("--output",metavar="DIR",help="Output folder the final matrix") 

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    outfile = open(str(args.output),"w")

    with open(str(args.input),"r") as pre_matrix:
        for lines in pre_matrix:
            aux=lines.split("\t")
            count=1
            #print(aux)
            for a in aux:
                if(count==1):
                    count-=-1
                    outfile.write(str(a)+"\t")
                    continue
                if(count%2==0):
                    if("\n" in a):
                        outfile.write(str(a))
                    else:
                        outfile.write(str(a)+"\t")
                count-=-1

if __name__ == "__main__":
    main()
