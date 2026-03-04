import os
import sys
import argparse
import shutil
import subprocess

def main ():

    parser = argparse.ArgumentParser(
        description="Maps Regprecise motif database agains the PRR fasta file",
        usage="python %(prog)s ")


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
        if os.path.exists(args.output_file):
            os.remove(args.output_file)
        
        motif_db_path = (args.motif_db).rstrip("/")        

        for line in motif_file:
            line = line.strip()
            mast_map = motif_db_path+"/"+line+".txt"
            cmd_job = ["mast", mast_map, args.prr_fasta, "-norc", "-nostatus", "-hit_list", "-best"]
            jobs.append(cmd_job)
            print(cmd_job)  # Print the job command (optional, to mimic STDERR print in Perl)

    # Execute each job command
    with open(args.output_file, "a") as outfile:  # Create the output file
        for job in jobs:
            subprocess.run(job, stdout=outfile, check=True)

if __name__ == "__main__":
    main()

