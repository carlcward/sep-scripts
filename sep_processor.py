import sys, os, subprocess, argparse, re
from collections import OrderedDict, namedtuple
import blast, ape_tools
from math import ceil
from copy import deepcopy
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
	try:
		return trans_dict[codon.upper()].split("|")[0]
	except:
		return "X"

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


def get_aa_seq(frame, start, stop, dna):
	fdna = dna[frame:]
	seq = translate_codons_to_string(get_codons(fdna))
	return seq[start:stop]

# flip the coordinates around the middle of the transcript if the match is on the minus strand
def pivot_coords(coords, transcript_start, transcript_end):
	new_coords = []
	length = transcript_end - transcript_start
	for coord in coords:
		new_coords.insert(0,[transcript_start + transcript_end - coord[1] + 1, transcript_start + transcript_end - coord[0]])
	return new_coords

# get a list of genomic coords for that are where the sep is actually transcribed from
def get_sep_trans_coords(transcript_coord,transcript,sep_start,sep_end,exon_data):
	out_string = ''
	#print sep_start, sep_end, len(transcript)
	regex = re.search("chr([0-9XY]{1,2}):([0-9]*)\-([0-9]*) strand=(.)",transcript_coord)

	if regex != None:
		regex = regex.groups()
		chromosome = regex[0]
		transcript_start = int(regex[1])
		transcript_end = int(regex[2])
		strand = regex[3]
		transcript_length = len(transcript)
		exon_data = filter(lambda e: e.strip() != '', exon_data.split(";"))
		exons = []
		for exon in exon_data:
			d = exon.split(",")
			exons.append( [int(d[0]),int(d[1])] )
		print exons
		start_index = 0
		stop_index = 0
		pos = 1
		#print exons 

		if strand == '-':
			nc = pivot_coords([[sep_start,sep_end]], 1, transcript_length)
			sep_start, sep_end = nc[0][0], nc[0][1]
		print "new" + str((sep_start,sep_end))

		for i in range(len(exons)):
			exon = exons[i]
			rel_start = pos
			rel_end = pos + (exon[1] - exon[0])
			#print exons[i], rel_start, rel_end
			print exons[i][0], exons[i][1], rel_start, rel_end
			if rel_start <= sep_start and rel_end >= sep_start:
				start_index = i
				exons[i][0] += sep_start - rel_start

			if rel_start <= sep_end and rel_end >= sep_end:
				stop_index = i
				exons[i][1] = exons[i][1] - (rel_end - sep_end)#exons[i][0] + (sep_end - pos) - 1
			
			pos = rel_end+1
		#print start_index,stop_index
		print start_index, stop_index
		sep_trans_coords = exons[start_index:stop_index+1]
		if strand == '-':
			print 'minus strand'
		#	sep_trans_coords = pivot_coords(sep_trans_coords, transcript_start, transcript_end)
		#print "TRANS",map(lambda c: (c[0] - coord_start,c[1] - coord_start, c[1]-c[0]) ,sep_trans_coords),sep_start,sep_end,sep_end-sep_start
		print sep_trans_coords
		total = 0
		for tc in sep_trans_coords:
			total += (tc[1]- tc[0])
			out_string += "chr%s:%s-%s;" % (chromosome, tc[0] , tc[1] )
		print total, total/3
		#if sep_trans_coords == []:
		#	out_string += transcript_coord.split(' ')[0]
	return out_string

# chekck see if there is a Kozak consensus start sequence in frame of the sep start, and before the stop
def check_kozak(start,stop,frame,dna):
	variants = set(['ATG','CTG','GTG','TTG','AAG','ACG','AGG','ATA','ATC','ATT'])
	fdna = dna[frame:]
	codons = get_codons(fdna)[start:]
	for i in range(stop-start):
		cdn = codons[i]
		if cdn in variants:
			if codons[i-1][0] == 'A' or codons[i-1][0] == 'G':
				if codons[i+1][0] == 'G':
					return i+start, cdn
	return None, ''

# returns all necessary information for a given sep
def process_sep(peptide,dna):
	koz_i = None
	# get the index and frame of the peptide
	index, frame = find_protein_in_dna(peptide,dna)
	if index == None:
		return (None, dna,(),(),(),(),(),(),'Protein not found, check anti-sense',(),(),(), ())
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
		return (None, None,(),(),(),(),(),(),'No downstream stop',(),(),(), ())
	else:
		# return: peptide, dna, frame, protein start, protein stop, marked sep start, marked sep stop, length, kozak index, kozak codon, kozak length
		if start_type == 'Non-AUG':
			koz_i, koz_cdn = check_kozak(start, index, frame, dna)
		aa_seq = get_aa_seq(frame, p_start, stop, dna)
		if koz_i == None:
			return (peptide, dna, frame, index, index + len(peptide), start, stop, stop-p_start, start_type, None, '', '', aa_seq)
		else:
			
			return (peptide, dna, frame, index, index + len(peptide), start, stop, stop-p_start, start_type, koz_i, koz_cdn, stop - koz_i, aa_seq)

# given sep data, and a list of sep data does the sep data "match" with all others in the list
# That is - has the same gene, location, start type, and length
def matching_sep(sep_data, sep_list):
	if len(sep_list) == 0:
		return 0
	for i in range(len(sep_list)):
		hit = sep_list[i]
		if hit.start_type == sep_data.start_type and hit.location == sep_data.location and hit.length == sep_data.length and hit.annotation.split("#")[0] == sep_data.annotation.split("#")[0]:
			return i
	return None


# create a dictionary of singly hit seps, put all questionably located ones into a another one
def put_unique_dict(udict, sdict, p_data):
	
	if 'Not Annotated' in p_data.annotation and len(filter(lambda pd: '#' in pd.annotation and 'not found' not in pd.annotation,udict[p_data.peptide] + sdict[p_data.peptide])) > 0:
		# Dont insert RNAseq data if there is already RefSeq Data there, (but do if it a special case (anti-sense etc))
		# RefSeq Data will necesseraily have to come in first as blast results are processed before normal
		pass
	else:
		# check if its in the suspect dict
		if len(sdict[p_data.peptide]) > 0:
			# and insert it into the correct location
			match_loc = matching_sep(p_data, sdict[p_data.peptide])
			if match_loc != None:
				sdict[p_data.peptide].insert(match_loc+1,DataTuple._make(('','',p_data.annotation,'','','',p_data.dna,'','','', '')))
			else:
				sdict[p_data.peptide].append(p_data)
		elif len(udict[p_data.peptide]) == 0:
			udict[p_data.peptide].append(p_data)
		else:
			# do the same for the merged dict, but if there becomes a conflict, move it to the sdict and delete it from the udict
			match_loc = matching_sep(p_data, udict[p_data.peptide])
			if match_loc == None:
				sdict[p_data.peptide] = udict[p_data.peptide] + [p_data]
				#print sdict[p_data.peptide], 'SDICT'
				udict[p_data.peptide] = []
			else:
				udict[p_data.peptide].insert(match_loc+1,DataTuple._make(('','',p_data.annotation,'','','',p_data.dna,'','','','')))
	
# write a results styled line to the output
def write_results_line(out_file_text,coord, peptide, annotation, location, start_type, sep_length, dna, sep_start, koz_cdn, koz_l, aa_seq):
	out_file = open(out_file_text, 'a')
	if len(annotation.split("#")) > 1:
		gene_id = annotation.split('#')[-1].strip()                                                
		annotation_url_string = "=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/nuccore/%s\"%s\"%s\")" % (gene_id,delim, ''.join(annotation.replace("#",",")))
		annotation = annotation_url_string
	
	out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (coord, peptide, annotation, location, start_type, sep_length, koz_cdn,koz_l,aa_seq, dna, sep_start))
	out_file.close()
def write_ape_file(out_path, file_name, dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop, koz_i, annotation):
	num = write_dict.setdefault(file_name,0) + 1
	if num == 1:
		num = ''
	
	if koz_i == None:
		koz_start, koz_end = 0,0
	else:
		koz_start = ape_tools.index_frame_to_loc(koz_i,frame)
		koz_end = koz_start+2
	write_dict[file_name] += 1
	out_file = open(out_path + file_name + str(num) + ".ape","w")
	
	out_file.write(ape_tools.create_simple_ape_map(dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop, koz_start, koz_end, annotation))
	out_file.close()
	


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Online Blast Search')
	parser.add_argument("sep_in", type=str)
	parser.add_argument("out", type=str)
	parser.add_argument("--blast", "-b", dest="blast", required=True, type=str,help='options: yes, no, both - both will use both RNAseq data and online data from BLAST hits')
	parser.add_argument('--officedelim','-o',dest='delim',action='store_true',help='use OpenOffice formula delimeters')
	args = parser.parse_args()


	in_file = open (args.sep_in, 'rU')
	out_path = args.out.rstrip('/') + "/"

	# clear the directory of ape files
	if os.path.isdir(out_path):
		print "Out directory exists"
	else:
		print "Out directory does not exist, Created!"
		subprocess.call(['mkdir', out_path])
	# the delimeter
	delim = ','
	if args.delim:
		delim = ';'
	# dictionary of results
	out_dict = OrderedDict()
	unique_out_dict = OrderedDict()
	suspect_out_dict = OrderedDict()
	DataTuple = namedtuple('DataTuple',['coord','peptide','annotation','location','start_type','length','dna','sep_trans_coords','koz_cdn','koz_len','aa_seq'])
	# collect all the peptide info into a list
	peptide_data = []
	for line in in_file:
		tup = tuple(line[:-1].split("\t"))
		if len(tup) == 4:
			tup = tup + ('',)
		peptide_data.append(tup)
	# DATA: coords, PCNUM, peptide, transcript
	peptide_data = filter(lambda p: p[2].strip() != '',peptide_data)
	queries = list(set(d[2] for d in peptide_data))

	# 
	blast_dict = {}
	
	unique_file_text = out_path + "results_merged.csv"
	open(unique_file_text, 'w')
	suspect_file_text = out_path + "results_multiple.csv"
	open(suspect_file_text, 'w')

	write_results_line(unique_file_text,'Coordinates','Peptide','Annotation','Location','Start Type','Length','RNASeq Transcript / RefSeq Trascript from BLAST','Sep Transcript Coords', 'Kozak Codon', 'Possible Length','Protein Sequence')
	write_results_line(suspect_file_text,'Coordinates','Peptide','Annotation','Location','Start Type','Length','RNASeq Transcript / RefSeq Trascript from BLAST','Sep Transcript Coords', 'Kozak Codon', 'Possible Length','Protein Sequence')
	
	if args.blast == 'yes':
		out_file_text = out_path + "results_blast.csv"
		xfile = blast.execute_blast_program(queries,'tblastn','refseq_rna')
		blast_dict = blast.parse_blast_results(xfile, 0, True)
		
		# only iterate through data with blast
		peptide_data = filter( lambda pd: pd[2] in blast_dict.keys(), peptide_data )
	elif args.blast == 'no':
		out_file_text = out_path + "results.csv"
		blast_dict = {}
	elif args.blast == 'both':
		out_file_text = out_path + "results.csv"
		xfile = blast.execute_blast_program(queries,'tblastn','refseq_rna')
		blast_dict = blast.parse_blast_results(xfile, 0, True)
	
	else:
		out_file_text = out_path + "results.csv"
		blast_dict = blast.parse_blast_results(open(args.blast,'r'), 0, True)
	
	for c,n,pep,d,e in peptide_data:
		out_dict[pep] = {
			'normal' : [],
			'blast' : []
		}
		unique_out_dict[pep] = []
		suspect_out_dict[pep] = []
	
	# write a little header and clear file
	open(out_file_text, 'w')
	write_results_line(out_file_text,'Coordinates','Peptide','Annotation','Location','Start Type','Length','RNASeq Transcript / RefSeq Trascript from BLAST','Sep Transcript Coords', 'Kozak Codon', 'Possible Length','Protein Sequence')



	for coord,pcnum,peptide,dna,exons in peptide_data:
		# lines are tab delimited: coord, pc_number, peptide, sequence
		cds_start, cds_stop = 0, 0
		
		
		# if the blast returned a hit
		
		
		# two options - the peptide had a blast hit or it didnt and if weve done it already
		if peptide in blast_dict.keys() and len(out_dict[peptide]['blast']) == 0:
			
			print "Processing: %s\t Blast Matches: %s" % (peptide, len(blast_dict[peptide]))
			for hit in blast_dict[peptide]:

				seq = hit[3]

				cds_range = hit[2].split("..")
				if len(cds_range) != 2:
					cds_start = 0
					cds_stop = 0
				else:
					cds_start = cds_range[0]
					cds_stop = cds_range[1]
				name = hit[0]
				file_name = name[name.rfind("(")+1:name.rfind(")")]
				annotation = file_name + "# " + hit[1].split("|")[3]
				
				# unfortunately we have to account for failed sep processing
				# due to either lack of downstream stop codon, or failure to find protein (anti-sense blast hit)
				# the start_type_or_message variable returns the start_type back or the failure message 
				out_peptide,blast_dna,frame,peptide_start,peptide_end,start,stop,sep_length,start_type_or_message, koz_i, koz_cdn, koz_l, aa_seq = process_sep(peptide, seq)
				if out_peptide:
					# get the location in the coding sequence, or show no coding sequence (also set cds_start, cds_stop)
					location,cds_start,cds_stop = ape_tools.calculate_location_in_protein(start, stop, frame, cds_start, cds_stop)
					# write the results file, with the RNAseq dna
					out_dict[peptide]['blast'].append(deepcopy((out_file_text, coord, peptide, annotation, location, start_type_or_message, sep_length, blast_dna, '',koz_cdn,koz_l,aa_seq)))
					dt = DataTuple._make((coord, peptide, annotation, location, start_type_or_message, sep_length, blast_dna, '',koz_cdn,koz_l, aa_seq))
					put_unique_dict(unique_out_dict, suspect_out_dict, dt)

					# write the ape file with the returned blast DNA
					
					write_ape_file(out_path, file_name, blast_dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop, koz_i, annotation)
				else:
					# if the search failed, write the results -- no need for ape map
					dt = DataTuple._make((coord, peptide, start_type_or_message + ", " + annotation,'', '','', blast_dna,'','','',''))
					# If it just a anti-sense, still put it in the dictionary, it has value
					if blast_dna != None:
						put_unique_dict(unique_out_dict,suspect_out_dict,dt)
					out_dict[peptide]['blast'].append(deepcopy((out_file_text, coord, peptide, start_type_or_message + ", " + annotation, '', '', '', blast_dna,'', '', '', '')))
		print "Processing: %s\t" % (peptide)
		# standard file and map creation without blast data
		if args.blast != 'yes':
			out_peptide,out_dna,frame,peptide_start,peptide_end,start,stop,sep_length,start_type_or_message, koz_i, koz_cdn, koz_l, aa_seq = process_sep(peptide, dna)
			annotation = 'Not Annotated'
			location = ''
			file_name = peptide[:4]
			cds_start = 0
			cds_stop = 0

			# check if it found a valid sep (i.e. has a downstream stop codon)
			if out_peptide:
				# use exons to get coordinates for the sep transcript
				sep_trans_coords = get_sep_trans_coords(coord,dna,ape_tools.index_frame_to_loc(start,frame),ape_tools.index_frame_to_loc(stop,frame),exons)
		
				out_dict[peptide]['normal'].append(deepcopy((out_file_text, coord, peptide, annotation, location, start_type_or_message, sep_length, dna, sep_trans_coords, koz_cdn, koz_l, aa_seq)))
				dt = DataTuple._make(( coord, peptide, annotation, location, start_type_or_message, sep_length, dna, sep_trans_coords, koz_cdn, koz_l, aa_seq))
				put_unique_dict(unique_out_dict, suspect_out_dict, dt)
				#write_results_line(out_file_text, coord, peptide, annotation, location, start_type, sep_length, dna, ape_tools.index_frame_to_loc(start,frame))
				write_ape_file(out_path, file_name, dna, frame, peptide_start, peptide_end, start, stop, sep_length, cds_start, cds_stop, koz_i, '')
			else:
				dt = DataTuple._make((coord, peptide, annotation, start_type_or_message + ", " + annotation, '', '','','','','',''))
				
				if out_dna != None:
					put_unique_dict(unique_out_dict,suspect_out_dict,dt)
				out_dict[peptide]['normal'].append(deepcopy((out_file_text, coord, peptide, start_type_or_message + ", " + annotation, '', '', '', dna,'', '', '','')))
	# create the file
	for peptide,data in out_dict.items():
		for line in data['normal']:
			write_results_line(*line)
		for line in data['blast']:
			write_results_line(*line)
	# create the merged file
	for peptide,data in unique_out_dict.items():
		for line in data:
			write_results_line(unique_file_text,*line)

	for peptide,data in suspect_out_dict.items():
		for line in data:
			write_results_line(suspect_file_text,*line)


	



