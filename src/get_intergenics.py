import sys
from Bio import SeqIO
from Bio.Seq import Seq

ptt_gene_file = open(sys.argv[1],"r")


intergenic_start = int(sys.argv[3])          #150
intergenic_stop = int(sys.argv[4])            #10)
default_intergenic_zone = int(sys.argv[5])   #50
minimal_intergenic_zone = default_intergenic_zone + intergenic_stop
gene_distance = int(sys.argv[6])+1             #50


def get_inter_left_start(current_tuple):
    value=int()
    if(int(current_tuple[0])>=minimal_intergenic_zone):
        value=int(current_tuple[0])-intergenic_start
        if(value<0):
            value=0
    return value

def get_inter_left_stop(current_tuple):
    return int(current_tuple[0])+intergenic_stop


def get_inter_right_start(current_tuple,keys):
    value2=int()
    if((contigs_size[keys]-(int(current_tuple[1])-int(minimal_intergenic_zone)))>=minimal_intergenic_zone):
        value2=int(current_tuple[1])+intergenic_start
        if(value2>=contigs_size[keys]):
            value2=contigs_size[keys]
    return value2
        

def get_inter_right_stop(current_tuple):
    return int(current_tuple[1])-intergenic_stop


#############################################################################################################

def get_paired_inter_right_start(current_tuple,keys,id_tuples):
    value2=int()
    index=contigs_with_gene_positions[keys].index(current_tuple)
    #print("the real tuple to calculate ",index,current_tuple)

    
    if(id_tuples==len(contigs_with_gene_positions[keys])-1):
        if((contigs_size[keys]-(int(current_tuple[1])-int(minimal_intergenic_zone)))>=minimal_intergenic_zone):
            value2=int(current_tuple[1])+intergenic_start
            if(value2>=contigs_size[keys]):
                value2=contigs_size[keys]
                #print("value1")
    if(id_tuples!=len(contigs_with_gene_positions[keys])-1):
        value2=int(current_tuple[1])+intergenic_start
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
            value2=int(current_tuple[1])+intergenic_start
            if(value2>=contigs_size[keys]):
                value2=contigs_size[keys]
    if(id_tuples!=len(new_minus_array_positions)-1):
        value2=int(current_tuple[1])+intergenic_start
        if(index==len(new_minus_array_positions)-1 and value2>=contigs_size[keys]):
            value2=contigs_size[keys]
        if(index<len(new_minus_array_positions)-1 and value2>int(new_minus_array_positions[index+1][0])):
            value2=int(new_minus_array_positions[index+1][0])-1
    return value2


def get_paired_inter_right_stop(current_tuple):
    return int(current_tuple[1])-intergenic_stop


def get_paired_inter_left_start(current_tuple,position):
    value2=int()
    #print("current_tuple ",current_tuple)
    index=contigs_with_gene_positions[keys].index(current_tuple)
    if(position==0):
        value2=int(current_tuple[0])-intergenic_start
        if(value2<0):
            value2=0
    if(position>0):
        value2=int(current_tuple[0])-intergenic_start
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
        value2=int(current_tuple[0])-intergenic_start
        if(value2<0):
            value2=0
    if(position>0):
        value2=int(current_tuple[0])-intergenic_start
        if(index>0 and value2<int(new_plus_array_positions[index-1][1])):
            value2=int(new_plus_array_positions[index-1][1])+1
        if(index==0 and value2<0):
            value2=0
    return value2



def get_paired_inter_left_stop(current_tuple):
    return int(current_tuple[0])+intergenic_stop





contigs_size=dict()
contigs_dna=dict()

with open(sys.argv[2]) as handle:
    for record in SeqIO.parse(handle,"fasta"):
        contigs_size[record.id]=len(record.seq)
        contigs_dna[record.id]=record.seq

#print(contigs_dna)
#print(contigs_dna["TARA_X000000368_G_scaffold272518_2"][643:1687])
#print(contigs_dna["TARA_X000000368_G_scaffold272518_2"][643:1687].reverse_complement())

#exit()

flag=False

contigs_with_gene_names = dict()
contigs_with_gene_positions = dict()
contigs_with_gene_strain = dict()

line_counter=1
for lines in ptt_gene_file:
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
            contigs_with_gene_strain.setdefault(aux[0],[])
            
        if(aux[0] in contigs_with_gene_positions):
            #print(lines)
            position=aux[1].split("..")
            position=list(position)
            position[0]=int(position[0])+1
            contigs_with_gene_positions[aux[0]].append(tuple(position))
            contigs_with_gene_names[aux[0]].append(aux[4])
            contigs_with_gene_strain[aux[0]].append(aux[2])
        #print(lines.replace("\n",""))
    line_counter+=1

global_operon_counter=0

for keys in contigs_with_gene_strain.keys():
##    print (keys)


##############   MULTI SIDE OPERON REPLICATION ##############



    if(contigs_with_gene_strain[keys].count("+")>=1 and contigs_with_gene_strain[keys].count("-")>=1):
        #print(contigs_with_gene_strain[keys])
        #print(contigs_with_gene_positions[keys])
        #print(contigs_with_gene_names[keys])
        secondary_operon_counter=0
        new_plus_array_positions=[]
        new_plus_array_names=[]
        new_plus_array_strain=[]
        
        new_minus_array_positions=[]
        new_minus_array_names=[]
        new_minus_array_strain=[]
      
    
        for ids, strain_value in enumerate(contigs_with_gene_strain[keys]):
            if(strain_value=="+"):
                new_plus_array_positions.append(contigs_with_gene_positions[keys][ids])
                new_plus_array_names.append(contigs_with_gene_names[keys][ids])
                new_plus_array_strain.append("+")
            if(strain_value=="-"):
                new_minus_array_positions.append(contigs_with_gene_positions[keys][ids])
                new_minus_array_names.append(contigs_with_gene_names[keys][ids])
                new_minus_array_strain.append("-")
        #print("ORIGINAL + ",new_plus_array_positions,new_plus_array_names)
        #print("ORIGINAL - ",new_minus_array_positions,new_minus_array_names)

### CASE 1 , SINGLE GENE + STRAIN , MULTITHREAD ###
        if(len(new_plus_array_strain)==1):
            if(new_plus_array_strain[0]=="+"):
                #print("OUTLAYERS == ",0,new_plus_array_positions[0][0])
                if((int(new_plus_array_positions[0][0]) + intergenic_stop ) >= minimal_intergenic_zone):
                    #print("only 1 + ",new_plus_array_positions,new_plus_array_names)
                    current_tuple=new_plus_array_positions[0]
                    inter_left_start=get_inter_left_start(current_tuple)
                    inter_left_stop=get_inter_left_stop(current_tuple)
                    contig_base_name=new_plus_array_names[0].split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strain[0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(new_plus_array_names[0]) )
                    #print("only 1 + ",str(new_plus_array_positions[0])+" side "+str(new_plus_array_strain[0])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop))
                    secondary_operon_counter+=1

### CASE 2 , SINGLE GENE - STRAIN , MULTITHREAD ###

        
        if(len(new_minus_array_strain)==1):
            if(new_minus_array_strain[0]=="-"):
                #print("OUTLAYERS == ",new_minus_array_positions[0][1],contigs_size[keys])
                if( contigs_size[keys] - (int(new_minus_array_positions[0][1])  - intergenic_stop) >= minimal_intergenic_zone):
                    #print("only 1 - ",new_minus_array_positions,new_minus_array_names)
                    current_tuple=new_minus_array_positions[0]
                    inter_right_start=get_inter_right_start(current_tuple,keys)
                    inter_right_stop=get_inter_right_stop(current_tuple)
                    contig_base_name=new_minus_array_names[0].split("_")
                    contig_base_name="_".join(contig_base_name[:-1])
                    print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strain[0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(new_minus_array_names[0]) )
                    #print("only 1 - ",str(new_minus_array_positions[0])+" side "+str(new_minus_array_strain[0])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start))
                    secondary_operon_counter+=1
### CASE 3 , MULTI GENE ALL + STRAIN , MULTITHREAD ###


        new_operon=[]
        position_calculation=0
        flag=False
        no_start=False
        if(len(new_plus_array_strain)>1):
            # this segments checks for all individual genes and a tuple containing their start and stop position, stored as a tuple,
            # we start a counter of the tuples at 0 and move on to the next one indicating the direcction of the straind + wich is left to right
            # the new_plus_array_* variables poses all the + tuples, no matter if they are part of the same operon or not. 
            if(new_plus_array_strain.count("+")==len(new_plus_array_strain)):
                #print(new_plus_array_strain, new_plus_array_strain.count("+"))
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
                    #print(current_tuple,current_tuple[0])
                    if (id_tuples+1<=(len(new_plus_array_strain)-1) and (int(new_plus_array_positions[id_tuples+1][0]) - int(current_tuple[1])) >=gene_distance):
                        #print(" -- OPERON BREAK -- ")
                        if(no_start==True):
                            no_start=False
                            id_tuples+=1
                            position_calculation=id_tuples
                            new_operon=[]
                            continue
                        #print("VALUE INNIT ",( int(new_plus_array_positions[0][0]) + intergenic_stop),  minimal_intergenic_zone)
                        if(id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                            #print(" *** can't have shit in los santos ",contigs_with_gene_positions[keys][0],contigs_size[keys])
                            position_calculation+=1
                            new_operon=[]
                            flag=True
                        if not (id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                            if(len(new_operon)>1):
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strain[0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                                #print(" / current_tuple and position_calculation # ", current_tuple,position_calculation)
                                #print(" / "+str(current_tuple)+" side "+str(new_plus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                                secondary_operon_counter+=1
                                new_operon=[]
                                no_start=False
                                position_calculation=id_tuples+1
                                #print(flag)
                            else:
                                #print(" // current_tuple and position_calculation # ", current_tuple,position_calculation)
                                inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                                inter_left_stop=get_paired_inter_left_stop(current_tuple)
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strain[0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                                #print(" // "+str(current_tuple)+" side "+str(new_plus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                                secondary_operon_counter+=1
                                new_operon=[]
                                position_calculation=id_tuples+1
                                no_start=False
                        if(id_tuples+1==(len(new_plus_array_strain)-1)):
                            flag=True


                    if ( id_tuples+1 <= (len(new_plus_array_strain)-1)  and  (int(new_plus_array_positions[id_tuples+1][0]) - int(current_tuple[1]) )<gene_distance):
                        #print(" --- contig calculation ---")
                        if(id_tuples < len(new_plus_array_strain)-1 and no_start==True):
                            id_tuples+=1
                            #print(" == CONTINUE ==")
                            continue

                        if(id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                            #print(" /// current_tuple and position_calculation - ", current_tuple,position_calculation)
                            #print(" /// can't have shit in los santos ",new_plus_array_positions[id_tuples],contigs_size[keys])
                            position_calculation+=1
                            new_operon=[]
                            no_start=True
                        if not (id_tuples == 0 and ( int(new_plus_array_positions[0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                            #print(" //// current_tuple and position_calculation  ", current_tuple,position_calculation)
                            current_tuple=new_plus_array_positions[position_calculation]
                            inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                            inter_left_stop=get_paired_inter_left_stop(current_tuple)
                            flag=False


                    if(id_tuples == len(new_plus_array_strain)-1):
                        #print(" -- END OF ARRAY --")
                        #print("no_start  =",no_start)
                        #print(" operon   =",new_operon)
                        #print("id_tuples =",id_tuples)
                        #print("position_calculation = ",position_calculation)
                        if(no_start==False and len(new_operon)==1 and flag==False):
                            flag=True

                        if(flag):
                            inter_left_start=get_paired_inter_left_start_multi(current_tuple,position_calculation,new_plus_array_positions)
                            inter_left_stop=get_paired_inter_left_stop(current_tuple)
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strain[0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            #print(" + end "+str(current_tuple)+" side "+str(new_plus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                        else:
                            if(len(new_operon)!=0 and no_start==False):
                                #print(" something something ")
                                #print(" - end "+str(current_tuple)+" side "+str(new_plus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_plus_array_strain[0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                                secondary_operon_counter+=1
                                1+1
                            if (len(new_operon)==1 and no_start==True):
                                #print(" # NO MAIDENS # ")
                                1+1
                            if (len(new_operon)==0 and no_start==True):
                                #print(" ## NO MAIDENS ## ")
                                1+1
                        new_operon=[]
                    id_tuples+=1
                #print("#### + END + #### ")



### CASE 4 , MULTI GENE ALL - STRAIN , MULTITHREAD ###


        new_operon=[]
        position_calculation=len(new_minus_array_strain)-1
        flag=False
        no_start=False
        if(len(new_minus_array_strain)>1):
            if(new_minus_array_strain.count("-")==len(new_minus_array_strain)):
                #print(new_minus_array_strain, new_minus_array_strain.count("-"))
                #print("###### ---------------- #############")
                #print(new_minus_array_positions)
                #print(new_minus_array_names)
                #print(" -- HARD LIMIT -- ", contigs_size[keys])
                #for id_tuples in range(len(new_minus_array_positions)-1,-1,-1):
                id_tuples=len(new_minus_array_positions)-1
                contig_base_name="".join(new_minus_array_names[0]).split("_")
                contig_base_name="_".join(contig_base_name[:-1])
                #secondary_operon_counter=0
                while (id_tuples >= 0):
                    #print("id_tuples and pos ",id_tuples,position_calculation)
                    #print(new_minus_array_positions[id_tuples])
                    current_tuple=new_minus_array_positions[id_tuples]
                    #print(new_minus_array_names)
                    
                    new_operon.insert(0,new_minus_array_names[id_tuples])

                    if (id_tuples-1>=0 and  (int(current_tuple[0])-int(new_minus_array_positions[id_tuples-1][1]))>=gene_distance):
                       #print("VALUE INNIT ",( (contigs_size[keys] - int(new_minus_array_positions[len(new_minus_array_strain)-1][1]) ) + intergenic_stop),  minimal_intergenic_zone)
                        if(no_start==True):
                            no_start=False
                            id_tuples-=1
                            position_calculation=id_tuples
                            new_operon=[]
                            continue
                        if(id_tuples == len(new_minus_array_positions)-1 and  contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strain)-1][1])  - intergenic_stop) < minimal_intergenic_zone):
                            #print("can't have shit in los santos ",new_minus_array_positions[len(new_minus_array_strain)-1],contigs_size[keys])
                            position_calculation-=1
                            new_operon=[]
                            flag=True
                        if not (id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strain)-1][1])  - intergenic_stop) < minimal_intergenic_zone):
                            if(len(new_operon)>1):
                                #print(" / "+str(current_tuple)+" side "+str(new_minus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_start)+" "+str(inter_right_stop)+" "+str(new_operon))
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strain[0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                                secondary_operon_counter+=1
                                new_operon=[]
                                no_start=False
                                position_calculation=id_tuples-1
                            else:
                                inter_right_start=get_paired_inter_right_start_multi(current_tuple,keys,position_calculation,new_minus_array_positions)
                                inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                #print(" // "+str(current_tuple)+" side "+str(new_minus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strain[0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                                secondary_operon_counter+=1
                                new_operon=[]
                                position_calculation=id_tuples-1
                                no_start=False
                        if(id_tuples-1==0):
                                flag=True

                    if ( no_start==False and id_tuples-1 >= 0 and (int(current_tuple[0])-int(new_minus_array_positions[id_tuples-1][1]))<gene_distance):
                        if(id_tuples > 0 and no_start==True):
                            id_tuples-=1
                            continue
                        if( id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strain)-1][1])  - intergenic_stop) < minimal_intergenic_zone):
                            #print("can't have shit in los santos ",new_minus_array_positions[len(new_minus_array_strain)-1],contigs_size[keys])
                            position_calculation-=1
                            new_operon=[]
                            no_start=True
                        if not ( id_tuples == len(new_minus_array_positions)-1 and contigs_size[keys] - (int(new_minus_array_positions[len(new_minus_array_strain)-1][1])  - intergenic_stop)  < minimal_intergenic_zone):
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
                            #print(" + end "+str(current_tuple)+" side "+str(new_minus_array_positions[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strain[0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                        
                        else:
                            if(len(new_operon)!=0 and no_start==False):
                                1+1
                                #print(" - end "+str(current_tuple)+" side "+str(new_minus_array_strain[position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(new_minus_array_strain[0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                                secondary_operon_counter+=1
                            if(len(new_operon)==1 and no_start==True):
                                #print(" # NO MAIDENS #")
                                1+1
                            if(len(new_operon)==0 and no_start==True):
                                #print(" ## NO MAIDENS ##")
                                1+1
                        new_operon=[]
                    id_tuples-=1

       

 #print(" -- END -- ")


##############   UNIQ SIDE OPERON REPLICATION ##############



### CASE 1 , SINGLE GENE + STRAIN ###


    if(len(contigs_with_gene_strain[keys])==1):
        #print(contigs_with_gene_positions[keys])
        if(contigs_with_gene_strain[keys][0]=="+"):
            if not ((int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                current_tuple=contigs_with_gene_positions[keys][0]
                inter_left_start=get_inter_left_start(current_tuple)
                inter_left_stop=get_inter_left_stop(current_tuple)
                contig_base_name=contigs_with_gene_names[keys][0].split("_")
                contig_base_name="_".join(contig_base_name[:-1])
                #secondary_operon_counter=0
                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_0,"+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(contigs_with_gene_names[keys][0]) )
                #print("only 1 - ",str(new_minus_array_positions[0])+" side "+str(new_minus_array_strain[0])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start))
                global_operon_counter+=1

                #print(str(contigs_with_gene_positions[keys][0])+" side "+str(contigs_with_gene_strain[keys][0])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop))


### CASE 2 , SINGLE GENE - STRAIN ###


    if(len(contigs_with_gene_strain[keys])==1):
        #print(contigs_with_gene_names[keys])
        if(contigs_with_gene_strain[keys][0]=="-"):
            if not ( contigs_size[keys] - (int(contigs_with_gene_positions[keys][0][1]) - intergenic_stop) < minimal_intergenic_zone):
                current_tuple=contigs_with_gene_positions[keys][0]
                inter_right_start=get_inter_right_start(current_tuple,keys)
                inter_right_stop=get_inter_right_stop(current_tuple)
                contig_base_name=contigs_with_gene_names[keys][0].split("_")
                contig_base_name="_".join(contig_base_name[:-1])
                #secondary_operon_counter=0
                print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_0,"+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(contigs_with_gene_names[keys][0]) )
                #print("only 1 - ",str(new_minus_array_positions[0])+" side "+str(new_minus_array_strain[0])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start))
                global_operon_counter+=1
                #print(str(contigs_with_gene_positions[keys][0])+" side "+str(contigs_with_gene_strain[keys][0])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start))


### CASE 3 , MULTI GENE ALL + STRAIN ###

    new_operon=[]
    position_calculation=0
    flag=False
    no_start=False
    if(len(contigs_with_gene_strain[keys])>1):
        if(contigs_with_gene_strain[keys].count("+")==len(contigs_with_gene_strain[keys])):
            #print(contigs_with_gene_strain[keys], contigs_with_gene_strain[keys].count("+"))
            #print("#######################################")
            #print(contigs_with_gene_positions[keys])
            #print(contigs_with_gene_names[keys])
            #for id_tuples in range(len(contigs_with_gene_positions[keys])):
            id_tuples=0
            contig_base_name="".join(contigs_with_gene_names[keys][0]).split("_")
            contig_base_name="_".join(contig_base_name[:-1])
                #secondary_operon_counter=0
            secondary_operon_counter=0
            while id_tuples < (len(contigs_with_gene_positions[keys])):
                # print(contigs_with_gene_positions[keys][id_tuples])
                #print(current_tuple,current_tuple[0])
                #print("id_tuples and pos ",id_tuples, position_calculation)
                new_operon.append(contigs_with_gene_names[keys][id_tuples])
                current_tuple=contigs_with_gene_positions[keys][id_tuples]
                if (id_tuples+1<=(len(contigs_with_gene_strain[keys])-1) and  (int(contigs_with_gene_positions[keys][id_tuples+1][0]) - int(current_tuple[1]))>=gene_distance):
                    #print(" -- OPERON BREAK -- ")
                    if(no_start==True):
                        no_start=False
                        id_tuples+=1
                        position_calculation=id_tuples
                        new_operon=[]
                        continue
                    #print("VALUE INNIT ",( int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop),  minimal_intergenic_zone)
                    if(id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                        #print(" *** can't have shit in los santos ",contigs_with_gene_positions[keys][0],contigs_size[keys])
                        position_calculation+=1
                        new_operon=[]
                        flag=True
                    if not (id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                        if(len(new_operon)>1):
                            #print(" / current_tuple and position_calculation # ", current_tuple,position_calculation)
                            #print(" / "+str(current_tuple)+" side "+str(contigs_with_gene_positions[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            new_operon=[]
                            no_start=False
                            position_calculation=id_tuples+1
                            #print(flag)
                        else:
                            #print(" // current_tuple and position_calculation # ", current_tuple,position_calculation)
                            inter_left_start=get_paired_inter_left_start(current_tuple,position_calculation)
                            inter_left_stop=get_paired_inter_left_stop(current_tuple)
                            #print(" // "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            new_operon=[]
                            position_calculation=id_tuples+1
                            no_start=False
                    if(id_tuples+1==(len(contigs_with_gene_strain[keys])-1)):
                        flag=True

                    
                if ( id_tuples+1 <= (len(contigs_with_gene_strain[keys])-1)  and  (int(contigs_with_gene_positions[keys][id_tuples+1][0]) - int(current_tuple[1]))<gene_distance):
                    #print(" --- contig calculation ---")
                    if(id_tuples < len(contigs_with_gene_strain[keys])-1 and no_start==True):
                        id_tuples+=1
                        #print(" == CONTINUE ==")
                        continue

                    if(id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                        #print(" /// current_tuple and position_calculation - ", current_tuple,position_calculation)
                        #print(" /// can't have shit in los santos ",contigs_with_gene_positions[keys][id_tuples],contigs_size[keys])
                        position_calculation+=1
                        new_operon=[]
                        no_start=True
                    if not (id_tuples == 0 and ( int(contigs_with_gene_positions[keys][0][0]) + intergenic_stop ) < minimal_intergenic_zone):
                        #print(" //// current_tuple and position_calculation  ", current_tuple,position_calculation)
                        current_tuple=contigs_with_gene_positions[keys][position_calculation]
                        inter_left_start=get_paired_inter_left_start(current_tuple,position_calculation)
                        inter_left_stop=get_paired_inter_left_stop(current_tuple)
                        flag=False
                    

                if(id_tuples == len(contigs_with_gene_strain[keys])-1):
                    #print(" -- END OF ARRAY --")
                    #print("no_start  =",no_start)
                    #print(" operon   =",new_operon)
                    #print("id_tuples =",id_tuples)
                    #print("flag =",flag)
                    #print("current_tuple =",current_tuple)
                    if(no_start==False and len(new_operon)==1 and flag==False):
                        flag=True
                    
                    if(flag):
                        inter_left_start=get_paired_inter_left_start(current_tuple,position_calculation)
                        inter_left_stop=get_paired_inter_left_stop(current_tuple)
                        #print(" + end "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                        print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                        secondary_operon_counter+=1
                    else:
                        if(len(new_operon)!=0 and no_start==False):
                            #print(" something something ")
                            #print(" - end "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_left_start)+" "+str(inter_left_stop)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_left_start:inter_left_stop]+","+str(inter_left_start)+","+str(inter_left_stop)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            1+1
                        if (len(new_operon)==1 and no_start==True):
                            #print(" # NO MAIDENS # ")
                            1+1
                        if (len(new_operon)==0 and no_start==True):
                            #print(" ## NO MAIDENS ## ")
                            1+1
                    new_operon=[]
                id_tuples+=1        


### CASE 4 , MULTI GENE ALL - STRAIN ###

    new_operon=[]
    position_calculation=len(contigs_with_gene_strain[keys])-1
    flag=False
    no_start=False
    
    if(len(contigs_with_gene_strain[keys])>1):
        if(contigs_with_gene_strain[keys].count("-")==len(contigs_with_gene_strain[keys])):
            
            #print("#######################################")
            #print(contigs_with_gene_strain[keys], contigs_with_gene_strain[keys].count("-"))
            #print(contigs_with_gene_positions[keys])
            #print(contigs_with_gene_names[keys])
            #print(" HARD LIMIT ",contigs_size[keys])
            #for id_tuples in range(len(contigs_with_gene_positions[keys])-1,-1,-1):
            id_tuples=len(contigs_with_gene_positions[keys])-1
            contig_base_name="".join(contigs_with_gene_names[keys][0]).split("_")
            contig_base_name="_".join(contig_base_name[:-1])
                #secondary_operon_counter=0
            secondary_operon_counter=0
            while (id_tuples >= 0):
                #print("id_tuples and pos ",id_tuples,position_calculation)
                #print(contigs_with_gene_positions[keys][id_tuples])
                current_tuple=contigs_with_gene_positions[keys][id_tuples]
                new_operon.insert(0,contigs_with_gene_names[keys][id_tuples])

                if (id_tuples-1>=0 and (int(current_tuple[0])-int(contigs_with_gene_positions[keys][id_tuples-1][1]))>=gene_distance):
                    #print("VALUE INNIT ",( (contigs_size[keys] - int(contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1][1]) ) + intergenic_stop),  minimal_intergenic_zone)
                    #print(" -- BREAK OPERON --")
                    if(no_start==True):
                        no_start=False
                        id_tuples-=1
                        position_calculation=id_tuples
                        new_operon=[]
                        continue
                    if(id_tuples == len(contigs_with_gene_positions[keys])-1 and  contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1][1]) - intergenic_stop) < minimal_intergenic_zone):
                        #print("can't have shit in los santos ",contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1],contigs_size[keys])
                        position_calculation-=1
                        new_operon=[]
                        flag=True
                    if not (id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1][1]) - intergenic_stop) < minimal_intergenic_zone):
                        if(len(new_operon)>1):
                            #print(" / "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_start)+" "+str(inter_right_stop)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            new_operon=[]
                            no_start=False
                            position_calculation=id_tuples-1
                        else:
                            inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                            inter_right_stop=get_paired_inter_right_stop(current_tuple)
                            #print(" // "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                            new_operon=[]
                            position_calculation=id_tuples-1
                            no_start=False

                    if(id_tuples-1==0):
                        flag=True



                if ( no_start==False and id_tuples-1 >= 0   and   (int(current_tuple[0])-int(contigs_with_gene_positions[keys][id_tuples-1][1]))<gene_distance  ):
                    #print(" ++ START OPERON ++")
                    if(id_tuples > 0 and no_start==True):
                        id_tuples-=1
                        continue
                    if( id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1][1]) - intergenic_stop) < minimal_intergenic_zone):
                        #print("can't have shit in los santos ",contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1],contigs_size[keys])
                        position_calculation-=1
                        new_operon=[]
                        no_start=True
                    if not ( id_tuples == len(contigs_with_gene_positions[keys])-1 and contigs_size[keys] - (int(contigs_with_gene_positions[keys][len(contigs_with_gene_strain[keys])-1][1]) - intergenic_stop) < minimal_intergenic_zone):
                        current_tuple=contigs_with_gene_positions[keys][position_calculation]
                        inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                        inter_right_stop=get_paired_inter_right_stop(current_tuple)
                        flag=False

                if(id_tuples == 0):
                    #print("final operon ",new_operon)
                    #print("id_tuples    ",id_tuples)
                    #print("position_calculation ",position_calculation)
                    #print("no start ",no_start)
                    #print("flag ",flag)
                    if(no_start==False and len(new_operon)==1 and flag==False):
                        flag=True
                    if(flag):
                        inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                        inter_right_stop=get_paired_inter_right_stop(current_tuple)
                        #print(" + end "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                        print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                        secondary_operon_counter+=1
                    else:
                        if(len(new_operon)!=0 and no_start==False):
                            1+1
                            if(len(new_operon)==1):
                                inter_right_start=get_paired_inter_right_start(current_tuple,keys,position_calculation)
                                inter_right_stop=get_paired_inter_right_stop(current_tuple)
                                #print("inside this shit ")
                            #print(" - end "+str(current_tuple)+" side "+str(contigs_with_gene_strain[keys][position_calculation])+" len:"+str(contigs_size[keys])+" # "+str(inter_right_stop)+" "+str(inter_right_start)+" "+str(new_operon))
                            print(str(contig_base_name)+",O_"+str(global_operon_counter)+"_"+str(secondary_operon_counter)+","+str(contigs_with_gene_strain[keys][0])+","+contigs_dna[keys][inter_right_stop:inter_right_start].reverse_complement()+","+str(inter_right_stop)+","+str(inter_right_start)+","+str(",".join(new_operon)) )
                            secondary_operon_counter+=1
                        if(len(new_operon)==1 and no_start==True):
                            #print(" # NO MAIDENS #")
                            1+1
                        if(len(new_operon)==0 and no_start==True):
                            #print(" ## NO MAIDENS ##")
                            1+1
                    new_operon=[]
                id_tuples-=1
    global_operon_counter+=1

