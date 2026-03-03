import sys
import argparse
import subprocess


def execute_calculate(args):
    pass

def execute_profiling(args):
    pass

def execute_matrix(args):
    pass



def main():

    parser = argparser.ArgumentParser(
        description="",
        usage="")

    subparsers = execute_calculate.add_subparsers(dest="command",required=True)

    execute_calculate = susparsers.add_parser("calculate",help ="")
    execute_calculate.add_argument("--sample_id",metavar="",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    execute_calculate.add_argument("--cds_file",metavar="",action=check_extension_pos({'ptt','gff','gff3'}),required=True,help="Path to file containing the positions of all CDS along your
  genomic sample. (PTT or GFF format)")
    execute_calculate.add_argument("--assembly",metavar="",action=check_extension_fasta({'fasta','fa','fna'}),required=True,help="Path to fasta formated genome")
    execute_calculate.add_argument("--prr_start",metavar="",type=int,default=300,required=True,help="Distance upstream of first CDS from operon start to consider for your PRR")
    execute_calculate.add_argument("--prr_stop",metavar="",type=int,default=30,required=True,help="Distance downtream of first CDS from operon start to consider for your PRR")
    execute_calculate.add_argument("--cds_dist",metavar="",type=int,default=50,required=True,help="Maximum distance between 2 CDS to consider to be part of the same CDS")
    execute_calculate.add_argument("--offset",metavar="",type=int,default=50,required=True,help="Distance from the edge of contigs to consider as offset")
    





if __name__ == "__main__":
    main()
