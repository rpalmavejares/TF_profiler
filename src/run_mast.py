from pathlib import Path
import os
import sys
import argparse
import shutil
import subprocess

def main ():

    parser = argparse.ArgumentParser(
        description="Maps Regprecise motif database agains the PRR fasta file",
        usage="python %(prog)s --motifs_list <FILE> --motif_db <DIR> --sample_id <ID> --prr_fasta <FILE> --output_folder <DIR> ")


    parser.add_argument("--motif_list",metavar="FILE",required=True,type=Path,default=Path("data/Regprecise_TF_DB/motifs.list"),help="List of motifs files to map agains the PRR fasta file")
    parser.add_argument("--motif_db",metavar="DIR",required=True,type=Path,default=Path("data/Regprecise_TF_DB/"),help="Path to Regprecise motifs DB")
    parser.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser.add_argument("--prr_fasta",metavar="FILE",required=True,help="Path to file containing PRR sequences")
    parser.add_argument("--output_folder",metavar="DIR",required=True,help="Output folder for all intermediate files")
    
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
        if os.path.exists(args.sample_id+"_prr.map"):
            os.remove(args.sample_id+"_prr.map")
        
        motif_db_path = str(args.motif_db).rstrip("/")        

        for line in motif_file:
            line = line.strip()
            mast_map = motif_db_path+"/"+line+".txt"
            cmd_job = ["mast", mast_map, args.prr_fasta, "-norc", "-nostatus", "-hit_list", "-best"]
            jobs.append(cmd_job)
            #print("Command Call: "+" ".join(cmd_job))
    # Execute each job command
    output_path = Path(args.output_folder) / f"{args.sample_id}_prr.map"
    with open(output_path, "a") as outfile:  # Create the output file
        for job in jobs:
            try:
                subprocess.run(job, stdout=outfile, check=True)
            except FileNotFoundError:
                print("Error: The program MAST from the MEME Suite has not been found")
            except subprocess.CalledProcessError as e:
                print(f"Error: The command ran but exited with an error code: {e.returncode}")

    outfile.close()

if __name__ == "__main__":
    main()

