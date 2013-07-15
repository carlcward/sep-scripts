import blast, sys
from fasta_reader import Data_Reader

queries = list(Data_Reader(sys.argv[1]))
blast_dict = blast.parse_blast_results(blast.execute_blast_program(queries,'blastp','nr'), 1, False)
results = [pep for pep in queries if pep not in blast_dict.keys()]
out_file = open(sys.argv[2], 'w')

for pep in queries:
	if pep in results:
		out_file.write(pep + "\tNOT FOUND\n")
	else:
		out_file.write(pep + "\tfound\n")