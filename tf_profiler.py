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

    call_command_calculate = ["python", "get_intergenics.py",
    "--sample_id", str(args.sample_id),
    "--cds_pos",str(args.cds_pos),
    "--assembly",str(args.assembly),
    "--prr_start",str(args.prr_start),
    "--prr_stop",str(args.prr_stop),
    "--cds_dist",str(args.cds_dist),
    "--offset",str(args.offset),
    "--output_folder",str(args.output_folder)]
    print(" ".join(call_command_calculate))
   
    folder_path = Path(args.output_folder)    
    folder_path.mkdir(parents=True, exist_ok=True) 
    subprocess.run(call_command_calculate,check=True)

def execute_profiling(args):
    pass

def execute_matrix(args):
    pass



def main():

    parser = argparse.ArgumentParser(
        description="",
        usage="")

    subparsers = parser.add_subparsers(dest="command",required=True)

    parser_c = subparsers.add_parser("calculate",help ="")
    parser_c.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser_c.add_argument("--cds_pos",metavar="FILE",action=check_extension_pos({'ptt','gff','gff3'}),required=True, help="Path to file containing the positions of all CDS along your genomic sample. (PTT or GFF format)")
    parser_c.add_argument("--assembly",metavar="FILE",action=check_extension_fasta({'fasta','fa','fna'}),required=True, help="Path to fasta formated genome")
    parser_c.add_argument("--prr_start",metavar="INT",type=int,default=300,help="Distance upstream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--prr_stop",metavar="INT",type=int,default=30,help="Distance downtream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--cds_dist",metavar="INT",type=int,default=50,help="Maximum distance between 2 CDS to consider to be part of the same CDS")
    parser_c.add_argument("--offset",metavar="INT",type=int,default=50,help="Distance from the edge of contigs to consider as offset")
    parser_c.add_argument("--output_folder",metavar="DIR",help="Output Folder for all intermediate files")
    
    parser_c.set_defaults(func=execute_calculate)    

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    main()
