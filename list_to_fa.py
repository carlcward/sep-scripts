#converts a line delimited list of sequences to a fasta with the sequence being the description
import sys
import argparse
def list_to_fa_string(seqs):
	out =''
	for seq in seqs:
		out += '>' + seq + '\n'
		out += seq + '\n'
	return out



if __name__ == '__main__':
	parser = argparse.ArgumentParser("Takes a line deliminted list and converts it to FASTA")
	parser.add_argument("in_file",type=str)
	parser.add_argument("out_file",type=str)
	args = parser.parse_args()
	in_file = open(args.in_file, "r")
	out_file = open(args.out_file, "w")
	#delchars = ''.join(c for c in map(chr, range(256)) if not c.isalpha())
	for line in in_file:
		line = line.strip()
		if line != '':
			clean = filter(str.isalpha,line)#line.translate(None, delchars)
			out_file.write(">" + clean + "\n")
			out_file.write(clean + "\n")