import sys
from fasta_reader import Data_Reader
from fa_index import check_peptide, load_indexes

out_file = open(sys.argv[3], 'w')
fasta, loc_file, ii_dict, length_dict = load_indexes(sys.argv[2])

for pep in Data_Reader(sys.argv[1]):
	res = check_peptide(pep, fasta, loc_file, ii_dict, length_dict)
	if res:
		out_file.write(pep + "\t" + res[0] + "\n")
	else:
		out_file.write(pep + "\t" + "NOT FOUND" + "\n")