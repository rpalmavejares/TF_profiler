import sys
import argparse
import subprocess
import os
from pathlib import Path
import glob

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
    print("All calculations completed at "+str(calculate_folder))



def execute_profiling(args):

    folder_path = Path(args.output_folder)    
    folder_path.mkdir(parents=True, exist_ok=True)
    output_folder = str(args.output_folder).rstrip("/")    
    profiling_folder = Path(output_folder+"/profiling")
    profiling_folder.mkdir(parents=True, exist_ok=True)

    run_profile_tf_feature= ["python", str(Path("src/profile_TF-Features.py")),
    "--sample_id",str(args.sample_id),
    "--annotation",str(args.annotation),
    "--feature_list",str(args.feature_list),
    "--coverage",str(args.coverage),
    "--cov_mode",str(args.cov_mode),
    "--cutoff",str(args.cutoff),
    "--tf_list",str(args.tf_list),
    "--prr_results",str(args.prr_results),
    "--output_folder",str(profiling_folder)]
    print("Command Call: "+" ".join(run_profile_tf_feature))
    subprocess.run(run_profile_tf_feature,check=True)
    print("Done")


    run_profile_tf= ["python", str(Path("src/profile_TF.py")),
    "--sample_id",str(args.sample_id),
    "--coverage",str(args.coverage),
    "--cov_mode",str(args.cov_mode),
    "--cutoff",str(args.cutoff),
    "--tf_list",str(args.tf_list),
    "--prr_results",str(args.prr_results),
    "--output_folder",str(profiling_folder)]
    print("Command Call: "+" ".join(run_profile_tf))
    subprocess.run(run_profile_tf,check=True)
    print("Done")


    run_profile_feature= ["python", str(Path("src/profile_Features.py")),
    "--sample_id",str(args.sample_id),
    "--annotation",str(args.annotation),
    "--feature_list",str(args.feature_list),
    "--coverage",str(args.coverage),
    "--cov_mode",str(args.cov_mode),
    "--output_folder",str(profiling_folder)]
    print("Command Call: "+" ".join(run_profile_feature))
    subprocess.run(run_profile_feature,check=True)
    print("Done")

    print("All profiles completed at "+str(profiling_folder))



def execute_matrix(args):

    folder_path = Path(args.output_folder)    
    folder_path.mkdir(parents=True, exist_ok=True)
    output_folder = str(args.output_folder).rstrip("/")    
    profiling_folder = Path(output_folder+"/profiling")
    matrix_folder = Path(output_folder+"/matrix")
    matrix_folder.mkdir(parents=True, exist_ok=True)
    working_cohort = output_folder.split("/")
    working_cohort = working_cohort[-1]


    ### PASTE PROFILES FEATURES ###

    p_features_folder = Path(str(profiling_folder)+"/P_Features")
    p_tf_folder = Path(str(matrix_folder)+"/P_TF")

    p_features_files = glob.glob(str(p_features_folder / "*.tsv"))    
    
    if(len(p_features_files) <=1):
        print("You need at least 2 profiles in the working folders to procced")
        sys.exit(1)

    temporal_feature_matrix = matrix_folder / "tmp_feature_matrix.tsv"

    run_paste_profiles_feature = ["paste"] + p_features_files + ["-d","\t"]

    try:
        with open(temporal_feature_matrix, "w") as out_file_feature:
            subprocess.run(run_paste_profiles_feature, stdout=out_file_feature, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
 
    ### MATRIX PROFILES FEATURES ###
 
    

    run_matrix_feature= ["python", str(Path("src/create_matrix.py")),
    "--input",str(temporal_feature_matrix),
    "--output",str(matrix_folder)+"/Matrix_"+str(working_cohort)+"_P_Features.tsv"]

    print("Command Call: "+" ".join(run_matrix_feature))
    subprocess.run(run_matrix_feature,check=True)
    print("Done")

   

    #print("Command Call: "+" ".join(run_paste_profiles))
    #subprocess.run(run_paste_profiles,check=True)
    #print("Done")
    
    #run_create_matrix=[]
    #print("Command Call: "+" ".join(run_create_matrix))
    #subprocess.run(run_create_matrix,check=True)
    #print("Done")


def main():

    parser = argparse.ArgumentParser(
        description="TF Profiler is a tools that computes a Simplified Operon Model to search for Potential Regulatory Regions and calculate profile abundances with MetaG and MetaT",
        usage="python %(prog)s <option> ")

    subparsers = parser.add_subparsers(dest="command",required=True)

    #### Calculate Menu Options ####

    parser_c = subparsers.add_parser("calculate",help ="Computes a Simplified Operon Model to extract Potential Regulatory Regions (PRR) on genomic contigs")
    parser_c.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser_c.add_argument("--cds_pos",metavar="FILE",action=check_extension_pos({'ptt','gff','gff3'}),required=True, help="Path to file containing the positions of all CDS along your genomic sample. (PTT or GFF format)")
    parser_c.add_argument("--assembly",metavar="FILE",action=check_extension_fasta({'fasta','fa','fna'}),required=True, help="Path to fasta formated genome")
    parser_c.add_argument("--prr_start",metavar="INT",type=int,default=300,help="Distance upstream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--prr_stop",metavar="INT",type=int,default=30,help="Distance downtream of first CDS from operon start to consider for your PRR")
    parser_c.add_argument("--cds_dist",metavar="INT",type=int,default=50,help="Maximum distance between 2 CDS to consider to be part of the same CDS")
    parser_c.add_argument("--offset",metavar="INT",type=int,default=50,help="Distance from the edge of contigs to consider as offset")
    parser_c.add_argument("--motif_list",metavar="FILE",type=Path,default=Path("data/Regprecise_TF_DB/motifs.list"),help="List of motifs files to map agains the PRR fasta file")
    parser_c.add_argument("--motif_db",metavar="DIR",type=Path,default=Path("data/Regprecise_TF_DB/"),help="Path to Regprecise motifs DB")
    parser_c.add_argument("--output_folder",metavar="DIR",help="Output folder for all project intermediary and results file")
    
    parser_c.set_defaults(func=execute_calculate)    


    #### Profile Menu Options ####

    parser_p = subparsers.add_parser("profiling",help ="Combines the PRR model with the Annotation and Coverage of CDS to comput a Regulation Profile on a single Sample")
    parser_p.add_argument("--sample_id",metavar="ID",required=True,help="Sample ID or Name. It must be contained on your CDS, Contigs, Coverage and Annotation nomenclature")
    parser_p.add_argument("--annotation",metavar="FILE",required=True,help="Annotation file containing your cds names and a feature (gene_name, COG, KO, EC_Number, etc. Must be a 2 column tab-separated file")
    parser_p.add_argument("--feature_list",metavar="FILE",required=True,help="")
    parser_p.add_argument("--coverage",metavar="FILE",required=True,help="")
    parser_p.add_argument("--cov_mode",metavar="MODE",choices=['cds','contig'],default="cds",required=True,help="Format of your genomic coverage. It can be used as 'Coverage per Contigs' or 'Coverage per CDS' ")
    parser_p.add_argument("--cutoff",metavar="FLOAT",type=float,default=1e6,required=True,help="E-value cutoff for mapped motifs agains PRR")
    parser_p.add_argument("--tf_list",metavar="FILE",default=str(Path("data/TF_Regprecise_list.txt")),help="List of all TFs to analize in the pipeline")
    parser_p.add_argument("--prr_results",metavar="FILE",type=Path,default=None,help="PRR Results file from the Calculate step")
    parser_p.add_argument("--output_folder",metavar="DIR",required=True,help="Output folder for all project intermediary and results file")

    parser_p.set_defaults(func=execute_profiling)    

    
    #### Matrix Menu Options ####
    
    
    parser_m = subparsers.add_parser("matrix",help ="Merge the output profiles of many samples or a cohort and convert the results into an abundance matrix")
    parser_m.add_argument("--output_folder",metavar="DIR",required=True,help="Output folder for all project intermediary and results file")
    
    parser_m.set_defaults(func=execute_matrix)    


    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if (args.command == "profiling"):
        if (args.prr_results is None):
            output_folder = str(args.output_folder).rstrip("/")    
            calculate_folder = Path(output_folder+"/calculate")
            args.prr_results = str(calculate_folder)+"/"+args.sample_id+"_prr.results"

    args.func(args)
    



if __name__ == "__main__":
    main()
