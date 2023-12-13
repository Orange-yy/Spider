# Find all spider toxins in uniprot and output them in a fasta file.
# First of all we are only interested in spider genes. 
# The means that the OC line must contain "Araneae".
# Observe the CC lines in the entries. 
# It can be seen that there are several sections and 
# they start with "-!-" and capital letters describing the section. 
# Pay attention to the following 3 sections: FUNCTION, TISSUE SPECIFICITY and SIMILARITY. 

import sys
import time

time_start = time.time()

# Get file name and open file
# filename = input("What swisprot file should I open:")
try:
    infile = open('uniprot_sprot.dat', 'r')
except IOError as err:
    print("Can't open file, reason", str(err))
    sys.exit(1)

# open the write file
try:
    outfile = open('toxin.fsa', 'w')
except IOError as err:
    print("Can't open file, reason", str(err))
    sys.exit(1)

# initialize variables used in loop for extraction
sp_id = ''
OCflag = False
toxin_flag = False
SQflag = False
aminoseq = ''
sum_count = 0
check_strings = ('toxin', 'toxic', 'venom', 'Venom', 'Toxin')
exclude_strings = ('Non-toxic', 'non-toxic', 'no toxic', 'not toxic', 'Not toxic', 'no toxicity', 'No toxicity')

# read line
for line in infile:
    # read line sign
    identifier = line[:2]

    # extract ID from the line
    if identifier == 'ID':
        # extract ID from the line
        start = 2
        # find start of id
        while line[start] == ' ':
            start += 1
        stop = start + 1
        # find end of id
        while line[stop] != ' ':
            stop += 1
        sp_id = line[start:stop]
        OCflag = False
        toxin_flag = False

    # Read the DE
    if identifier == 'DE':
        if any(s in line for s in check_strings):
            toxin_flag = True
        if any(s in line for s in exclude_strings):
            toxin_flag = False

    # Read the OC 
    if identifier == 'OC':
        # check the species
        if 'Araneae' in line:
            OCflag = True
            # initialize variables
            function_count = 0
            function_flag = False
            tissue_count = 0
            tissue_flag = False
            similarity_count = 0
            similarity_flag = False

    # deal with CC lines
    if OCflag and identifier == 'CC':

        # deal with function lines
        if line[5:8] == '-!-' and function_flag:
            function_flag = False
        if line[5:17] == '-!- FUNCTION':
            function_flag = True
        if function_flag:
            if any(s in line for s in exclude_strings):
                function_count = -3
                function_flag = False
            elif any(s in line for s in check_strings):
                function_count = 1
                function_flag = False
            else:
                function_count = 0

        # deal with tissue specificity lines
        if line[5:8] == '-!-' and tissue_flag:
            tissue_flag = False
        if line[5:27] == '-!- TISSUE SPECIFICITY':
            tissue_flag = True
        if tissue_flag:
            if 'venom gland' in line:
                tissue_count = 1
                tissue_flag = False
            else:
                tissue_count = 0

        # deal with tissue similarity lines
        if line[5:8] == '-!-' and similarity_flag:
            similarity_flag = False
        if line[5:19] == '-!- SIMILARITY':
            similarity_flag = True
        if similarity_flag:
            if any(s in line for s in check_strings):
                similarity_count = 1
                similarity_flag = False
            else:
                similarity_count = 0

        # plus every count
        sum_count = function_count + tissue_count + similarity_count

        # bonus score
        if line[5:19] == '-!- TOXIC DOSE':
            sum_count += 1
        if toxin_flag:
            sum_count += 1

    # bonus score
    if OCflag and identifier == 'KW':
        if 'Toxin' in line:
            sum_count += 2

    # Getting sequence using stateful parsing
    # Red line
    if identifier == '//' and sum_count >= 3:
        SQflag = False
        sum_count = 0
    # Collect data
    if SQflag:
        aminoseq += line
    # Green line
    if identifier == 'SQ' and sum_count >= 3:
        SQflag = True

    # deal with SQ
    if aminoseq != '' and not SQflag:
        # Remove spaces, treatment of data after extraction.
        cleanseq = ''
        for i in range(len(aminoseq)):
            if aminoseq[i] != ' ':
                cleanseq = cleanseq + aminoseq[i]
        aminoseq = ''

        # Output
        outfile.write(">" + sp_id + '\n')
        # output amino seq
        for k in range(0, len(cleanseq), 60):
            outfile.write(cleanseq[k:k+60])


infile.close()
outfile.close()
time_end = time.time()
print('time cost', time_end-time_start, 's')


