import os
import sys

# Get command-line arguments
mast = sys.argv[1]
pathdb = sys.argv[2]
pathout = sys.argv[3]
motif = sys.argv[4]
fasta = sys.argv[5]

# Open the motif file for reading
with open(motif, 'r') as motif_file:
    jobs = []
    output = f"{pathout}.txt"
    
    # Remove the output file if it exists and create a new one
    if os.path.exists(output):
        os.remove(output)
    open(output, 'a').close()  # Create the output file

    for line in motif_file:
        line = line.strip()
        meme = os.path.join(pathdb, f"{line}.txt")
        stringjob = f"{mast} {meme} {fasta} -norc -nostatus -hit_list -best >> {output}"
        jobs.append(stringjob)
        print(stringjob)  # Print the job command (optional, to mimic STDERR print in Perl)

# Execute each job command
for job in jobs:
    os.system(job)
