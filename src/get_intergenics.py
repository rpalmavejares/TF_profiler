import sys
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
from pathlib import Path

def get_inter_left_start(current_tuple):
    value=int()
    if(int(current_tuple[0])>=minimal_intergenic_zone):
        value=int(current_tuple[0])-prr_start
        if(value<0):
            value=0
    return value

def get_inter_left_stop(current_tuple):
    return int(current_tuple[0])+prr_stop


def get_inter_right_start(current_tuple,keys):
    value2=int()
    if((contigs_size[keys]-(int(current_tuple[1])-int(minimal_intergenic_zone)))>=minimal_intergenic_zone):
        value2=int(current_tuple[1])+prr_start
        if(value2>=contigs_size[keys]):
            value2=contigs_size[keys]
    return value2
        

def get_inter_right_stop(current_tuple):
    return int(current_tuple[1])-prr_stop


#############################################################################################################

def get_paired_inter_right_start(current_tuple,keys,id_tuples):
    value2=int()
    index=contigs_with_gene_positions[keys].index(current_tuple)
    #print("the real tuple to calculate ",index,current_tuple)

    
    if(id_tuples==len(contigs_with_gene_positions[keys])-1):
        if((contigs_size[keys]-(int(current_tuple[1])-int(minimal_intergenic_zone)))>=minimal_intergenic_zone):
            value2=int(current_tuple[1])+prr_start
            if(value2>=contigs_size[keys]):
                value2=contigs_size[keys]
                #print("value1")
    if(id_tuples!=len(contigs_with_gene_positions[keys])-1):
        value2=int(current_tuple[1])+prr_start
        if(index==len(contigs_with_gene_positions[keys])-1 and value2>=contigs_size[keys]):
            value2=contigs_size[keys]
            #print("value2")
            return value2
        if(index<=len(contigs_with_gene_positions[keys])-1 and value2>int(contigs_with_gene_positions[keys][index+1][0])):
            value2=int(contigs_with_gene_positions[keys][index+1][0])-1
            #print("value3")
    return value2

def get_paired_inter_right_start_multi(current_tuple,keys,id_tuples,new_minus_array_positions):
    value2=int()
    index=new_minus_array_positions.index(current_tuple)
    #print("the real tuple to calculate ",index,current_tuple)


    if(id_tuples==len(new_minus_array_positions)-1):
        if((contigs_size[keys]-(int(current_tuple[1])-int(minimal_intergenic_zone)))>=minimal_intergenic_zone):
            value2=int(current_tuple[1])+prr_start
            if(value2>=contigs_size[keys]):
                value2=contigs_size[keys]
    if(id_tuples!=len(new_minus_array_positions)-1):
        value2=int(current_tuple[1])+prr_start
        if(index==len(new_minus_array_positions)-1 and value2>=contigs_size[keys]):
            value2=contigs_size[keys]
        if(index<len(new_minus_array_positions)-1 and value2>int(new_minus_array_positions[index+1][0])):
            value2=int(new_minus_array_positions[index+1][0])-1
    return value2


def get_paired_inter_right_stop(current_tuple):
    return int(current_tuple[1])-prr_stop


def get_paired_inter_left_start(keys,current_tuple,position):
    value2=int()
    #print("current_tuple ",current_tuple)
    index=contigs_with_gene_positions[keys].index(current_tuple)
    if(position==0):
        value2=int(current_tuple[0])-prr_start
        if(value2<0):
            value2=0
    if(position>0):
        value2=int(current_tuple[0])-prr_start
        if(index>0 and value2<int(contigs_with_gene_positions[keys][index-1][1])):
            value2=int(contigs_with_gene_positions[keys][index-1][1])+1
        if(index==0 and value2<0):
            value2=0            
    return value2

def get_paired_inter_left_start_multi(current_tuple,position,new_plus_array_positions):
    value2=int()
    #print("current_tuple ",current_tuple)
    index=new_plus_array_positions.index(current_tuple)
    if(position==0):
        value2=int(current_tuple[0])-prr_start
        if(value2<0):
            value2=0
    if(position>0):
        value2=int(current_tuple[0])-prr_start
        if(index>0 and value2<int(new_plus_array_positions[index-1][1])):
            value2=int(new_plus_array_positions[index-1][1])+1
        if(index==0 and value2<0):
            value2=0
    return value2



def get_paired_inter_left_stop(current_tuple):
    return int(current_tuple[0])+prr_stop


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

def main():


    parser = argparse.ArgumentParser(
        description="",
        usage="")

# Output == <Sample_Name>_operon_model.csv

    
    parser.add_argument("--sample_id",metavar="",required=True,help="Sample ID or Name. It must be contained on your CDS and Contigs nomenclature")
    parser.add_argument("--cds_pos",metavar="",action=check_extension_pos({'ptt','gff','gff3'}),required=True,help="Path to file containing the positions of all CDS along your genomic sample. (PTT or GFF format)")
    parser.add_argument("--assembly",metavar="",action=check_extension_fasta({'fasta','fa','fna'}),required=True,help="Path to fasta formated genome")
    parser.add_argument("--prr_start",metavar="",type=int,default=300,required=True,help="Distance upstream of first CDS from operon start to consider for your PRR")
    parser.add_argument("--prr_stop",metavar="",type=int,default=30,required=True,help="Distance downtream of first CDS from operon start to consider for your PRR")
    parser.add_argument("--cds_dist",metavar="",type=int,default=50,required=True,help="Maximum distance between 2 CDS to consider to be part of the same CDS")
    parser.add_argument("--offset",metavar="",type=int,default=50,required=True,help="Distance from the edge of contigs to consider as offset")

    if len(sys.argv) == 1:
        print("\n")
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    global cds_pos
    global prr_start, prr_stop, offset, minimal_intergenic_zone, cds_dist
    global contigs_size, contigs_dna 
    global contigs_with_gene_names
    global contigs_with_gene_positions
    global contigs_with_gene_strand
    
    prr_start = args.prr_start
    prr_stop = args.prr_stop
    offset = args.offset
    minimal_intergenic_zone = offset + prr_stop
    cds_dist = args.cds_dist + 1


    contigs_size=dict()
    contigs_dna=dict()

    with open(args.assembly) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            contigs_size[record.id]=len(record.seq)
            contigs_dna[record.id]=record.seq
    handle.close()

    flag=False

    contigs_with_gene_names = dict()
    contigs_with_gene_positions = dict()
    contigs_with_gene_strand = dict()

    line_counter=1
    with open(args.cds_pos) as cpf:
        for lines in cpf:
            aux=lines.split("\t")
            if(line_counter==2 or line_counter==3):
                line_counter+=1
                continue
            if("Marine metagenome" in lines or "length=" in lines):
                line_counter=1
            if(line_counter>=4):
                if(aux[0] not in contigs_with_gene_positions):
                    contigs_with_gene_positions.setdefault(aux[0],[])
                    contigs_with_gene_names.setdefault(aux[0],[])
                    contigs_with_gene_strand.setdefault(aux[0],[])
                    
                if(aux[0] in contigs_with_gene_positions):
                    position=aux[1].split("..")
                    position=list(position)
                    position[0]=int(position[0])+1
                    contigs_with_gene_positions[aux[0]].append(tuple(position))
                    contigs_with_gene_names[aux[0]].append(aux[4])
                    contigs_with_gene_strand[aux[0]].append(aux[2])
            line_counter+=1

        global_operon_counter=0
    cpf.close()

    out_prr_model = Path(args.sample_id+"_prr_model.csv")
    out_prr_model = open(args.sample_id+"_prr_model.csv", "w")

    for keys in contigs_with_gene_strand.keys():


    ##############   MULTI SIDE OPERON REPLICATION ##############



        if(contigs_with_gene_strand[keys].count("+")>=1 and contigs_with_gene_strand[keys].count("-")>=1):
            secondary_operon_counter=0
            new_plus_array_positions=[]
            new_plus_array_names=[]
            new_plus_array_strand=[]
            
            new_minus_array_positions=[]
            new_minus_array_names=[]
            new_minus_array_strand=[]
          
        
            for ids, strand_value in enumerate(contigs_with_gene_strand[keys]):
                if(strand_value=="+"):
                    new_plus_array_positions.append(contigs_with_gene_positions[keys][ids])
                    new_plus_array_names.append(contigs_with_gene_names[keys][ids])
                    new_plus_array_strand.append("+")
                if(strand_value=="-"):
                    new_minus_array_positions.append(contigs_with_gene_positions[keys][ids])
                    new_minus_array_names.append(contigs_with_gene_names[keys][ids])
                    new_minus_array_strand.append("-")


    ### CASE 1 , SINGLE GENE + STRAIN , MULTITHREAD ###


            if(len(new_plus_array_strand)==1):
                if(new_plus_array_strand[0]=="+"):
                    if((int(new_plus_array_positions[0][0]) + prr_stop ) >= minimal_intergenic_zone):
                        current_tuple=new_plus_array_positions[0]
                        inter_left_start=get_inter_left_start(current_tuple)
                        inter_left_stop=get_inter_left_stop(current_tuple)
                        contig_base_name=new_plus_array_names[0].split("_")
                        contig_base_name="_".join(contig_base_name[:-1])
                        out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strand[0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(new_plus_array_names[0]+"\n") )
                        secondary_operon_counter+=1

    ### CASE 2 , SINGLE GENE - STRAIN , MULTITHREAD ###

            
            if(len(new_minus_array_strand)==1):
                if(new_minus_array_strand[0]=="-"):
                    if( contigs_size[keys] - (int(new_minus_array_positions[0][1])  - prr_stop) >= minimal_intergenic_zone):
                        current_tuple=new_minus_array_positions[0]
                        inter_right_start=get_inter_right_start(current_tuple,keys)
                        inter_right_stop=get_inter_right_stop(current_tuple)
                        contig_base_name=new_minus_array_names[0].split("_")
                        contig_base_name="_".join(contig_base_name[:-1])
                        out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strand[0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(new_minus_array_names[0])+"\n")
                        secondary_operon_counter+=1
    ### CASE 3 , MULTI GENE ALL + STRAIN , MULTITHREAD ###


            new_operon=[]
            position_calculation=0
            flag=False
            no_start=False
            if(len(new_plus_array_strand)>1):
                # this segments checks for all individual genes and a tuple containing their start and stop position, stored as a tuple,
                # we start a counter of the tuples at 0 and move on to the next one indicating the direcction of the strandd + wich is left to right
                # the new_plus_array_* variables poses all the + tuples, no matter if they are part of the same operon or not. 
                if(new_plus_array_strand.count("+")==len(new_plus_array_strand)):
                    #print(new_plus_array_strand, new_plus_array_strand.count("+"))
                    #print("### +++++++++++++++++++++++++++++ ###")
                    #print(new_plus_array_positions)
                    #print(new_plus_array_names)
                    #for id_tuples in range(len(new_plus_array_positions)):
                    #    print(id_tuples)
                    id_tuples=0
                    contig_base_name="".join(new_plus_array_names[0]).split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    # we will extract the name of the parent contig that contains the gene, so we can extract the respective intergenic region sequence from there.

                    while id_tuples <= (len(new_plus_array_positions)-1):
                        #print(id_tuples)
                        #print(new_plus_array_positions)
                        #print(new_plus_array_positions[id_tuples])
                        #print("id_tuples and pos ",id_tuples, position_calculation)
                        
                        # while you still have tuples, you continue to build operons, this is the start of a new one therefore the new_operon array
                        
                        new_operon.append(new_plus_array_names[id_tuples])
                        current_tuple=new_plus_array_positions[id_tuples]
                        if (id_tuples+1<=(len(new_plus_array_strand)-1) and (int(new_plus_array_positions[id_tuples+1][0]) - int(current_tuple[1])) >=cds_dist):
                            if(no_start==True):
                                no_start=False
                                id_tuples+=1
                                position_calculation=id_tuples
                                new_operon=[]
                                continue
                            if(id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + prr_stop ) < minimal_intergenic_zone):
                                position_calculation+=1
                                new_operon=[]
                                flag=True
                            if not (id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + prr_stop ) < minimal_intergenic_zone):
                                if(len(new_operon)>1):
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strand[0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                    new_operon=[]
                                    no_start=False
                                    position_calculation=id_tuples+1
                                else:
                                    inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                                    inter_left_stop=get_paired_inter_left_stop(current_tuple)
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strand[0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                    new_operon=[]
                                    position_calculation=id_tuples+1
                                    no_start=False
                            if(id_tuples+1==(len(new_plus_array_strand)-1)):
                                flag=True


                        if ( id_tuples+1 <= (len(new_plus_array_strand)-1)  and  (int(new_plus_array_positions[id_tuples+1][0]) - int(current_tuple[1]) )<cds_dist):
                            if(id_tuples < len(new_plus_array_strand)-1 and no_start==True):
                                id_tuples+=1
                                continue

                            if(id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + prr_stop ) < minimal_intergenic_zone):
                                position_calculation+=1
                                new_operon=[]
                                no_start=True
                            if not (id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + prr_stop ) < minimal_intergenic_zone):
                                current_tuple=new_plus_array_positions[position_calculation]
                                inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                                inter_left_stop=get_paired_inter_left_stop(current_tuple)
                                flag=False


                        if(id_tuples == len(new_plus_array_strand)-1):
                            if(no_start==False and len(new_operon)==1 and flag==False):
                                flag=True
                            if(flag):
                                inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                                inter_left_stop=get_paired_inter_left_stop(current_tuple)
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strand[0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                            else:
                                if(len(new_operon)!=0 and no_start==False):
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strand[0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                if (len(new_operon)==1 and no_start==True):
                                    pass
                                if (len(new_operon)==0 and no_start==True):
                                    pass
                            new_operon=[]
                        id_tuples+=1
                    #print("#### + END + #### ")



    ### CASE 4 , MULTI GENE ALL - STRAIN , MULTITHREAD ###


            new_operon=[]
            position_calculation=len(new_minus_array_strand)-1
            flag=False
            no_start=False
            if(len(new_minus_array_strand)>1):
                if(new_minus_array_strand.count("-")==len(new_minus_array_strand)):
                    id_tuples=len(new_minus_array_positions)-1
                    contig_base_name="".join(new_minus_array_names[0]).split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    while (id_tuples >= 0):
                        current_tuple=new_minus_array_positions[id_tuples]
                        new_operon.insert(0,new_minus_array_names[id_tuples])
                        if (id_tuples-1>=0 and  (int(current_tuple[0])-int(new_minus_array_positions[id_tuples-1][1]))>=cds_dist):
                            if(no_start==True):
                                no_start=False
                                id_tuples-=1
                                position_calculation=id_tuples
                                new_operon=[]
                                continue
                            if(id_tuples == len(new_minus_array_positions)-1 and  contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strand)-1][1])  - prr_stop) < minimal_intergenic_zone):
                                position_calculation-=1
                                new_operon=[]
                                flag=True
                            if not (id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strand)-1][1])  - prr_stop) < minimal_intergenic_zone):
                                if(len(new_operon)>1):
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strand[0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                    new_operon=[]
                                    no_start=False
                                    position_calculation=id_tuples-1
                                else:
                                    inter_right_start=get_paired_inter_right_start_multi(current_tuple,keys,position_calculation,new_minus_array_positions)
                                    inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strand[0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                    new_operon=[]
                                    position_calculation=id_tuples-1
                                    no_start=False
                            if(id_tuples-1==0):
                                    flag=True

                        if ( no_start==False and id_tuples-1 >= 0 and (int(current_tuple[0])-int(new_minus_array_positions[id_tuples-1][1]))<cds_dist):
                            if(id_tuples > 0 and no_start==True):
                                id_tuples-=1
                                continue
                            if( id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strand)-1][1])  - prr_stop) < minimal_intergenic_zone):
                                position_calculation-=1
                                new_operon=[]
                                no_start=True
                            if not ( id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strand)-1][1])  - prr_stop)  < minimal_intergenic_zone):
                                current_tuple=new_minus_array_positions[position_calculation]
                                inter_right_start=get_paired_inter_right_start_multi(current_tuple,keys,position_calculation,new_minus_array_positions)
                                inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                flag=False

                        if(id_tuples == 0):
                            if(no_start==False and len(new_operon)==1 and flag==False):
                                flag=True

                            if(flag):
                                inter_right_start=get_paired_inter_right_start_multi(current_tuple,keys,position_calculation,new_minus_array_positions)
                                inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strand[0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                            
                            else:
                                if(len(new_operon)!=0 and no_start==False):
                                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strand[0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                    secondary_operon_counter+=1
                                if(len(new_operon)==1 and no_start==True):
                                    pass
                                if(len(new_operon)==0 and no_start==True):
                                    pass
                            new_operon=[]
                        id_tuples-=1




    ##############   UNIQ SIDE OPERON REPLICATION ##############



    ### CASE 1 , SINGLE GENE + STRAIN ###


        if(len(contigs_with_gene_strand[keys])==1):
            if(contigs_with_gene_strand[keys][0]=="+"):
                if not ((int(contigs_with_gene_positions[keys][0][0]) + prr_stop ) < minimal_intergenic_zone):
                    current_tuple=contigs_with_gene_positions[keys][0]
                    inter_left_start=get_inter_left_start(current_tuple)
                    inter_left_stop=get_inter_left_stop(current_tuple)
                    contig_base_name=contigs_with_gene_names[keys][0].split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_0,"+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(contigs_with_gene_names[keys][0])+"\n")
                    global_operon_counter+=1



    ### CASE 2 , SINGLE GENE - STRAIN ###


        if(len(contigs_with_gene_strand[keys])==1):
            if(contigs_with_gene_strand[keys][0]=="-"):
                if not ( contigs_size[keys] - (int(contigs_with_gene_positions[keys][0][1]) - prr_stop) < minimal_intergenic_zone):
                    current_tuple=contigs_with_gene_positions[keys][0]
                    inter_right_start=get_inter_right_start(current_tuple,keys)
                    inter_right_stop=get_inter_right_stop(current_tuple)
                    contig_base_name=contigs_with_gene_names[keys][0].split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_0,"+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(contigs_with_gene_names[keys][0])+"\n")
                    global_operon_counter+=1


    ### CASE 3 , MULTI GENE ALL + STRAIN ###

        new_operon=[]
        position_calculation=0
        flag=False
        no_start=False
        if(len(contigs_with_gene_strand[keys])>1):
            if(contigs_with_gene_strand[keys].count("+")==len(contigs_with_gene_strand[keys])):
                id_tuples=0
                contig_base_name="".join(contigs_with_gene_names[keys][0]).split("_")
                contig_base_name="_".join(contig_base_name[:-1])
                secondary_operon_counter=0
                while id_tuples < (len(contigs_with_gene_positions[keys])):
                    new_operon.append(contigs_with_gene_names[keys][id_tuples])
                    current_tuple=contigs_with_gene_positions[keys][id_tuples]
                    if (id_tuples+1<=(len(contigs_with_gene_strand[keys])-1) and  (int(contigs_with_gene_positions[keys][id_tuples+1][0]) - int(current_tuple[1]))>=cds_dist):
                        if(no_start==True):
                            no_start=False
                            id_tuples+=1
                            position_calculation=id_tuples
                            new_operon=[]
                            continue
                        if(id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + prr_stop ) < minimal_intergenic_zone):
                            position_calculation+=1
                            new_operon=[]
                            flag=True
                        if not (id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + prr_stop ) < minimal_intergenic_zone):
                            if(len(new_operon)>1):
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                                new_operon=[]
                                no_start=False
                                position_calculation=id_tuples+1
                            else:
                                inter_left_start=get_paired_inter_left_start(keys,current_tuple,position_calculation)
                                inter_left_stop=get_paired_inter_left_stop(current_tuple)
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                                new_operon=[]
                                position_calculation=id_tuples+1
                                no_start=False
                        if(id_tuples+1==(len(contigs_with_gene_strand[keys])-1)):
                            flag=True

                        
                    if ( id_tuples+1 <= (len(contigs_with_gene_strand[keys])-1)  and  (int(contigs_with_gene_positions[keys][id_tuples+1][0]) - int(current_tuple[1]))<cds_dist):
                        if(id_tuples < len(contigs_with_gene_strand[keys])-1 and no_start==True):
                            id_tuples+=1
                            continue

                        if(id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + prr_stop ) < minimal_intergenic_zone):
                            position_calculation+=1
                            new_operon=[]
                            no_start=True
                        if not (id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + prr_stop ) < minimal_intergenic_zone):
                            current_tuple=contigs_with_gene_positions[keys][position_calculation]
                            inter_left_start=get_paired_inter_left_start(keys,current_tuple,position_calculation)
                            inter_left_stop=get_paired_inter_left_stop(current_tuple)
                            flag=False
                        

                    if(id_tuples == len(contigs_with_gene_strand[keys])-1):
                        if(no_start==False and len(new_operon)==1 and flag==False):
                            flag=True
                        if(flag):
                            inter_left_start=get_paired_inter_left_start(keys,current_tuple,position_calculation)
                            inter_left_stop=get_paired_inter_left_stop(current_tuple)
                            out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                            secondary_operon_counter+=1
                        else:
                            if(len(new_operon)!=0 and no_start==False):
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_left_start:inter_left_stop])+","+str(inter_left_start+1)+","+str(inter_left_stop+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                                1+1
                            if (len(new_operon)==1 and no_start==True):
                                pass
                            if (len(new_operon)==0 and no_start==True):
                                pass
                        new_operon=[]
                    id_tuples+=1        


    ### CASE 4 , MULTI GENE ALL - STRAIN ###

        new_operon=[]
        position_calculation=len(contigs_with_gene_strand[keys])-1
        flag=False
        no_start=False
        
        if(len(contigs_with_gene_strand[keys])>1):
            if(contigs_with_gene_strand[keys].count("-")==len(contigs_with_gene_strand[keys])):
                id_tuples=len(contigs_with_gene_positions[keys])-1
                contig_base_name="".join(contigs_with_gene_names[keys][0]).split("_")
                contig_base_name="_".join(contig_base_name[:-1])
                secondary_operon_counter=0
                while (id_tuples >= 0):
                    current_tuple=contigs_with_gene_positions[keys][id_tuples]
                    new_operon.insert(0,contigs_with_gene_names[keys][id_tuples])

                    if (id_tuples-1>=0 and (int(current_tuple[0])-int(contigs_with_gene_positions[keys][id_tuples-1][1]))>=cds_dist):
                        if(no_start==True):
                            no_start=False
                            id_tuples-=1
                            position_calculation=id_tuples
                            new_operon=[]
                            continue
                        if(id_tuples == len(contigs_with_gene_positions[keys])-1 and  contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strand[keys])-1][1]) - prr_stop) < minimal_intergenic_zone):
                            position_calculation-=1
                            new_operon=[]
                            flag=True
                        if not (id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strand[keys])-1][1]) - prr_stop) < minimal_intergenic_zone):
                            if(len(new_operon)>1):
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                                new_operon=[]
                                no_start=False
                                position_calculation=id_tuples-1
                            else:
                                inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                                inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                                new_operon=[]
                                position_calculation=id_tuples-1
                                no_start=False

                        if(id_tuples-1==0):
                            flag=True



                    if ( no_start==False and id_tuples-1 >= 0   and   (int(current_tuple[0])-int(contigs_with_gene_positions[keys][id_tuples-1][1]))<cds_dist  ):
                        if(id_tuples > 0 and no_start==True):
                            id_tuples-=1
                            continue
                        if( id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strand[keys])-1][1]) - prr_stop) < minimal_intergenic_zone):
                            position_calculation-=1
                            new_operon=[]
                            no_start=True
                        if not ( id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strand[keys])-1][1]) - prr_stop) < minimal_intergenic_zone):
                            current_tuple=contigs_with_gene_positions[keys][position_calculation]
                            inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                            inter_right_stop=get_paired_inter_right_stop(current_tuple)
                            flag=False

                    if(id_tuples == 0):
                        if(no_start==False and len(new_operon)==1 and flag==False):
                            flag=True
                        if(flag):
                            inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                            inter_right_stop=get_paired_inter_right_stop(current_tuple)
                            out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                            secondary_operon_counter+=1
                        else:
                            if(len(new_operon)!=0 and no_start==False):
                                if(len(new_operon)==1):
                                    inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                                    inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                out_prr_model.write(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strand[keys][0])+","+str(contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement())+","+str(inter_right_stop+1)+","+str(inter_right_start+1)+","+str(",".join(new_operon))+"\n")
                                secondary_operon_counter+=1
                            if(len(new_operon)==1 and no_start==True):
                                pass
                            if(len(new_operon)==0 and no_start==True):
                                pass
                        new_operon=[]
                    id_tuples-=1
        global_operon_counter+=1
    out_prr_model.close()
if __name__ == "__main__":
    main()

