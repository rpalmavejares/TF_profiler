import os
import sys
import argparse
import shutil
import subprocess

def main ():

    parser = argparse.ArgumentParser(
        description="",
        usage="")


    parser.add_argument("--motif_list",metavar="",help="")
    parser.add_argument("--motif_db",metavar="",help="")
    parser.add_argument("--sample_id",metavar="",help="")
    parser.add_argument("--prr_fasta",metavar="",help="")
    parser.add_argument("--output_file",metavar="",help="")
    
    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()


    # Get command-line arguments

    # Open the motif file for reading
    with open(args.motif_list, 'r') as motif_file:
        jobs = []
        # Remove the output file if it exists and create a new one
        if os.path.exists(output):
            os.remove(output)
        open(args.output_file+".txt", 'a').close()  # Create the output file
        
        motif_db_path = (args.motif_db).rstrip("/")        

        for line in motif_file:
            line = line.strip()
            mast_map = motif_db_path+"/"+line+".txt"
            stringjob = f"mast mast_map {args.prr_fasta} -norc -nostatus -hit_list -best >> {args.output_file}.txt"
            jobs.append(stringjob)
            print(stringjob)  # Print the job command (optional, to mimic STDERR print in Perl)

    # Execute each job command
    for job in jobs:
        subprocess.run(job)

if __name__ == "__main__":
    main()

