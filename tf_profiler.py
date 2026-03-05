import sys
import argparse
import subprocess
import os
from pathlib import Path

def check_extension_fasta(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, fname, option_string=None):
            # Extract the extension
            ext = os.path.splitext(fname)[1][1:].lower()
            if ext not in choices:
                parser.error(f"File must be the correct format: {choices}")
            setattr(namespace, self.dest, fname)
    return Act


def check_extension_pos(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, fname, option_string=None):
            # Extract the extension
            ext = os.path.splitext(fname)[1][1:].lower()
            if ext not in choices:
                parser.error(f"File must be the correct format: {choices}")
            setattr(namespace, self.dest, fname)
    return Act


def execute_calculate(args):


    folder_path = Path(args.output_folder)    
    folder_path.mkdir(parents=True, exist_ok=True)
    output_folder = str(args.output_folder).rstrip("/")    
    calculate_folder = Path(output_folder+"/calculate")
    calculate_folder.mkdir(parents=True, exist_ok=True)
    call_command_calculate_a = ["python", str(Path("src/get_intergenics.py")),
    "--sample_id", str(args.sample_id),
    "--cds_pos",str(args.cds_pos),
    "--assembly",str(args.assembly),
    "--prr_start",str(args.prr_start),
    "--prr_stop",str(args.prr_stop),
    "--cds_dist",str(args.cds_dist),
    "--offset",str(args.offset),
    "--output_folder",str(calculate_folder)]
    print("############################")
    print("### Running TF Profiler ####")
    print("############################")
    print("Command Call: "+" ".join(call_command_calculate_a))
    #print(" ".join(call_command_calculate))
    subprocess.run(call_command_calculate_a,check=True)
    print("Done")
    out_operon_model = (str(calculate_folder)+"/"+args.sample_id+"_operon_model.csv")   
    call_command_calculate_b = ["python", str(Path("src/get_fasta.py")),
    "--sample_id",str(args.sample_id),
    "--operon_model",str(out_operon_model),
    "--output_folder",str(calculate_folder)]
    print("Command Call: "+" ".join(call_command_calculate_b))
    subprocess.run(call_command_calculate_b,check=True)
    print("Done")
    output_prr_fasta = str(calculate_folder)+"/"+args.sample_id+"_prr.fasta"
    call_command_calculate_c = ["python", str(Path("src/run_mast.py")),
    "--motif_list",str(args.motif_list),
    "--motif_db",str(args.motif_db),
    "--sample_id",str(args.sample_id),
    "--prr_fasta",str(output_prr_fasta),
    "--output_folder",str(calculate_folder)]
    print("Command Call: "+" ".join(call_command_calculate_c))
    subprocess.run(call_command_calculate_c,check=True)
    print("Done")
    output_prr_map = str(calculate_folder)+"/"+args.sample_id+"_prr.map"
    output_prr_results = str(calculate_folder)+"/"+args.sample_id+"_prr.results"
    call_command_calculate_d = ["python", str(Path("src/parse_results.py")),
    "--sample_id",str(args.sample_id),
    "--prr_fasta",str(output_prr_fasta),
    "--motif_alignment",str(output_prr_map),
    "--output",str(output_prr_results)]
    print("Command Call: "+" ".join(call_command_calculate_c))
    subprocess.run(call_command_calculate_d,check=True)
    print("Done")



def execute_profiling(args):

    folder_path = Path(args.output_folder)    
    folder_path.mkdir(parents=True, exist_ok=True)
    output_folder = str(args.output_folder).rstrip("/")    
    profiling_folder = Path(output_folder+"/profiling")
    profiling_folder.mkdir(parents=True, exist_ok=True)






def execute_matrix(args):
    pass



def main():

    parser = argparse.ArgumentParser(
        description="",
        usage="")

    subparsers = parser.add_subparsers(dest="command",required=True)

    #### Calculate Menu Options ####

    parser_c = subparsers.add_parser("calculate",help ="")
    parser_c.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser_c.add_argument("--cds_pos",metavar="FILE",action=check_extension_pos({'ptt','gff','gff3'}),required=True, help="Path to file containing the positions of all CDS along your genomic sample. (PTT or GFF format)")
    parser_c.add_argument("--assembly",metavar="FILE",action=check_extension_fasta({'fasta','fa','fna'}),required=True, help="Path to fasta formated genome")
    parser_c.add_argument("--prr_start",metavar="INT",type=int,default=300,help="Distance upstream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--prr_stop",metavar="INT",type=int,default=30,help="Distance downtream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--cds_dist",metavar="INT",type=int,default=50,help="Maximum distance between 2 CDS to consider to be part of the same CDS")
    parser_c.add_argument("--offset",metavar="INT",type=int,default=50,help="Distance from the edge of contigs to consider as offset")
    parser_c.add_argument("--motif_list",metavar="FILE",type=Path,default=Path("data/Regprecise_TF_DB/motifs.list"),help="List of motifs files to map agains the PRR fasta file")
    parser_c.add_argument("--motif_db",metavar="DIR",type=Path,default=Path("data/Regprecise_TF_DB/"),help="Path to Regprecise motifs DB")
    parser_c.add_argument("--output_folder",metavar="DIR",help="Output folder for all intermediate files")
    
    parser_c.set_defaults(func=execute_calculate)    


    #### Calculate Menu Options ####


    parser_p = subparsers.add_parser("profiling",help ="")
    parser_p.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")


    parser_p.set_defaults(func=execute_profiling)    


    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    main()
