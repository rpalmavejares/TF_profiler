import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



intergenic_region_file = open(sys.argv[1],"r")

all_records = []

for lines in intergenic_region_file:
    genes=list()
    aux = lines.replace("\n","")
    aux=aux.split(",")
    replicon = aux[0]
    operon_name = aux[1]
    strand = aux[2]
    sequence = aux[3]
    start = aux[4]
    stop = aux[5]
    genes = aux[6:]
    if(len(sequence)!=0):
        description="Replicon="+str(replicon)+" Strand="+str(strand)+" Start="+str(start)+" Stop="+str(stop)+" Genes="+str(",".join(genes))

        #print(operon_name)
        record = SeqRecord(Seq(sequence), id=operon_name, description=description)
        all_records.append(record)

SeqIO.write(all_records,sys.argv[2],"fasta")


