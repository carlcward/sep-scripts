import sys, os, subprocess, argparse
import blast, ape_tools

# keep track of file numbers
write_dict = {}
# codon -> protein dict
trans_dict = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser", 
"TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp", 
"TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu", 
"CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro", 
"CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg", 
"CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met", 
"ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn", 
"AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg", 
"GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala", 
"GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu", 
"GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}

# get a list of codons given a sequence
def get_codons(seq):
	codons = []
	while len(seq) > 2:
		cdn = seq [:3]
		codons.append(cdn)
		seq = seq[3:]
	return codons

def reverse_compliment(seq):
	comp = ''
	base_comp = {'A':'T','C':'G','G':'C','T':'A'}
	for base in seq:
		comp += base_comp[base]
	# return reverse compliment
	return comp[::-1]

# translate a single codon to a single amino acid
def translate_codon(codon):
	return trans_dict[codon.upper()].split("|")[0]

# translate list of codons into list of 1 letter amino acid codes
def translate_codons(codons):
	return reduce (lambda a,b: a + [translate_codon(b)], codons, [] )

# convert list of codons to protein string
def translate_codons_to_string(codons):
	return ''.join(translate_codons(codons))

# NOT USED, translate given dna to one letter amino acid protein string
def translate_dna_seq_to_string(seq):
	return reduce (lambda a,b: a + translate_codon(b), get_codons(seq), "")

# get the location of the protein, as well as the frame (0,1,2) to shift the dna by
def find_protein_in_dna(protein, dna):
	frames = []
	frames.append((dna,0))
	frames.append((dna[1:],1))
	frames.append((dna[2:],2))
	for seq,frame in frames:
		pt = translate_codons_to_string(get_codons(seq))
		index = pt.find(protein)
		if index > -1:
			#return DNA locations
			return (index, frame)
	
	return (None, None)

# start at the end of the protein and find the first stop
def find_downstream_stop(index, frame, dna):
	fdna = dna[frame:]
	icdns = get_codons(fdna)
	stops = ['TAA','TGA','TAG']

	for i in range(index,len(icdns)):# index; i < len(icdns); i++:
		if icdns[i] in stops:
			return i 

	return None

# find the upstream starting features
def find_upstream_start(index, frame, dna):
	fdna = dna[frame:]
	# the indexed codons
	icdns = get_codons(fdna)
	stops = ['TAA','TGA','TAG']
	closest_stop = None
	furthest_start = None
	# from the start of the sequence until the protein starts
	for i in range(0, index):
		# if there's a stop, mark it as the closest and reset the start
		if icdns[i] in stops:
			closest_stop = i
			furthest_start = None
		# if there's a start, mark it at the farthest unless there already is one
		if furthest_start == None and icdns[i] == 'ATG':
			furthest_start = i
	return (furthest_start, closest_stop)

# returns all necessary information for a given sep
def process_sep(peptide,dna):
	# get the index and frame of the peptide
	index, frame = find_protein_in_dna(peptide,dna)
	if index == None:
		return (None, 'Protein not found, check anti-sense',(),(),(),(),(),(),())
	# find nearest downstream stop
	stop = find_downstream_stop(index,frame,dna)
	# find farthest upstream start, without a stop between
	start, u_stop = find_upstream_start(index,frame,dna)

	# where to start COUNTING for the peptide (upstream stops should NOT factor into length)
	p_start = start
	
	# determine case
	if start == None:
		start_type = 'Non-AUG'
		if u_stop == None:
			p_start = 0
			start = 0
		else:
			start = u_stop
			p_start = u_stop + 1
	else:
		start_type = 'AUG'
	if stop == None:
		return (None, 'No downstream stop',(),(),(),(),(),(),())
	else:
		# return: peptide, dna, frame, protein start, protein stop, marked sep start, marked sep stop, length
		return (peptide, dna, frame, index, index + len(peptide), start, stop, stop-p_start, start_type)


# write a results styled line to the output
def write_results_line(out_file,coord, peptide, annotation, location, start_type, sep_length, dna, sep_start):
	if len(annotation.split(",")) > 1:
		gene_id = annotation.split(',')[-1].strip()                                                
		annotation_url_string = "=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/nuccore/%s\";\"%s\")" % (gene_id, annotation)
		annotation = annotation_url_string
	out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (coord, peptide, annotation, location, start_type, sep_length, dna, sep_start))

def write_ape_file(out_path, file_name, dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop):
	num = write_dict.setdefault(file_name,0) + 1
	if num == 1:
		num = ''
	write_dict[file_name] += 1
	out_file = open(out_path + file_name + str(num) + ".ape","w")
	out_file.write(ape_tools.create_ape_map(dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop))
	out_file.close()
	


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Online Blast Search')
	parser.add_argument("sep_in", type=str)
	parser.add_argument("out", type=str)
	parser.add_argument("--blast", "-b", dest="blast", required=True, type=str)
	args = parser.parse_args()


	in_file = open (args.sep_in, 'rU')
	out_path = args.out.rstrip('/') + "/"
	# clear the directory of ape files
	if os.path.isdir(out_path):
		print "Out directory exists"
	else:
		print "Out directory does not exist, Created!"
		subprocess.call(['mkdir', out_path])
	
	
	# collect all the peptide info into a list
	peptide_data = []
	for line in in_file:
		peptide_data.append(tuple(line.strip().split("\t")))

	# DATE: coords, PCNUM, peptide, transcript
	peptide_data = filter(lambda p: p[0] != '',peptide_data)
	
	queries = set(d[2] for d in peptide_data)
	'''
	for d in peptide_data:
		if d[2] not in queries:
			queries.append(d[2])
	print queries
	'''
	# tblastn the peptides and collect the data
	if args.blast == 'yes':
		out_file_text = open(out_path + "results_blast.csv", 'w')
		req = blast.send_blast_request(queries)
		xml_data = blast.get_blast_results(req)
		blast_dict = blast.parse_blast_results(xml_data)
		# only iterate through data with blast
		peptide_data = filter( lambda pd: pd[2] in blast_dict.keys(), peptide_data )
	elif args.blast == 'no':
		out_file_text = open(out_path + "results.csv", 'w')
		blast_dict = {}
	else:
		out_file_text = open(out_path + "results_blast.csv", 'w')
		blast_dict = blast.parse_blast_results(open(args.blast,'r').read())
		peptide_data = filter( lambda pd: pd[2] in blast_dict.keys(), peptide_data )
	# write a little header
	write_results_line(out_file_text,'Coordinates','Peptide','Annotation','Location','Start Type','Length','RNASeq Transcript','Sep Start')
	
	blast_done = []

	for coord,pcnum,peptide,dna in peptide_data:
		# lines are tab delimited: coord, pc_number, peptide, sequence
		cds_start, cds_stop = 0, 0
		annotation = 'Not Annotated'
		
		file_name = peptide[:4]
		# if the blast returned a hit
		
		location = ''
		# two options - the peptide had a blast hit or it didnt
		# we have slightly different processing/iterator patterns, so we have to control for this unfortunately 
		if peptide in blast_dict.keys() :
			# check if we have already processed the blast data for this peptide
			if peptide not in blast_done:
				print "Processing: %s\t Blast Matches: %s" % (peptide, len(blast_dict[peptide]))
				blast_done.append(peptide)			
				for hit in blast_dict[peptide]:

					seq = hit[3]
					cds_range = hit[2].split("..")
					cds_start = cds_range[0]
					cds_stop = cds_range[1]
					name = hit[0]
					file_name = name[name.rfind("(")+1:name.rfind(")")]
					annotation = file_name + ", " + hit[1].split("|")[3]
					
					# unfortunately we have to account for failed sep processing
					# due to either lack of downstream stop codon, or failure to find protein (anti-sense blast hit)
					# the blast_dna_or_message variable returns the dna back (for the map) or the failure message 
					out_peptide,blast_dna_or_message,frame,peptide_start,peptide_end,start,stop,sep_length,start_type = process_sep(peptide, seq)
					
					if out_peptide:
						# get the location in the coding sequence, or show no coding sequence (also set cds_start, cds_stop)
						location,cds_start,cds_stop = ape_tools.calculate_location_in_protein(start, stop, cds_start, cds_stop, ape_tools.index_frame_to_loc(start,frame))
						# write the results file, with the RNAseq dna
						write_results_line(out_file_text, coord, peptide, annotation, location, start_type, sep_length, dna)
						# write the ape file with the returned blast DNA
						write_ape_file(out_path, file_name, blast_dna_or_message, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop)
					else:
						# if the search failed, write the results -- no need for ape map
						write_results_line(out_file_text, coord, peptide, blast_dna_or_message + ", " + annotation, '', '', '', dna,'')
		else:
			print "Processing: %s\t" % (peptide)
			# standard file and map creation without blast data
			out_peptide,out_dna_or_message,frame,peptide_start,peptide_end,start,stop,sep_length,start_type = process_sep(peptide, dna)
			# check if it found a valid sep (i.e. has a downstream stop codon)
			if out_peptide:
				write_results_line(out_file_text, coord, peptide, annotation, location, start_type, sep_length, dna, ape_tools.index_frame_to_loc(start,frame))
				write_ape_file(out_path, file_name, out_dna_or_message, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop)
			else:
				write_results_line(out_file_text, coord, peptide, out_dna_or_message, '', '', '', dna, '')
		
		

	out_file_text.close()



