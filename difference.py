import sys

#  open file
try:
    standardfile = open('uniprot-taxonomy%3Aaraneae+keyword%3Atoxin+AND+reviewed%3Ayes.fasta', 'r')
    myfile = open('toxin.fsa', 'r')
except IOError as err:
    print("Can't open reason", str(err))
    sys.exit(1)

# initializing varibles
standardset = set()
myset = set()
count1 = 0
count2 = 0

# get standardfile line
for line in standardfile:
    if line.startswith('>'):
        line = line.split()
        name = line[0]
        name = name.split('|')
        name = name[-1]
        standardset.add(name)
        count1 += 1
standardfile.close()

# get myfile line
for line in myfile:
    if line.startswith('>'):
        line = line.strip()
        line = line.strip('>')
        myset.add(line)
        count2 += 1
myfile.close()

# find difference
# diffset = standardset.difference(resset)
diffset = myset.difference(standardset)
difflist = list(diffset)
difflist.sort()

# output result
print("There are", count1, "protein in standardfile")
print("There are", count2, "protein in myfile")
print("myfile have", len(diffset), "different proteins")
print("They are")
print(difflist)