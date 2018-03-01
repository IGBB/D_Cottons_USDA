import glob
from Bio.Nexus import Nexus

# adapted from http://biopython.org/wiki/Concatenate_nexus
# the combine function takes a list of tuples [(name, nexus instance)...],
# if we provide the file names in a list we can use a list comprehension to
# create these tuples

file_list = glob.glob('*.nxs')
nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]

combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('all.nex', 'w'))