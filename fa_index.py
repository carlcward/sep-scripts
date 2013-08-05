import cPickle, re
def get_locations(iloc, loc_file):
	char ='g'
	data_string = ''
	loc_file.seek(iloc)
	while char != "#" and char != '':
		char = loc_file.read(1)
		data_string += char
	data_string =data_string[:-1]
	data_string = data_string.split(",") 
	return data_string[1:]

def check_loc(sloc,peptide,fasta,length_dict):
	fasta.seek(int(sloc))
	seq = fasta.read(length_dict[int(sloc)])
	if peptide in seq:
		fasta.seek(int(sloc))
		desc,char = '',''
		i = 0
		while char != '>':
			i += 1
			fasta.seek(int(sloc)-i)
			char = fasta.read(1)
			if char != '\n':
				desc += char
		return desc[::-1],seq
	return None

# essentially a boolean -> was the peptide found?
def check_peptide(peptide, fasta, loc_file, ii_dict, length_dict):
	key = peptide[:5]
	if key in ii_dict:
		iloc = int(ii_dict[peptide[:5]])
	else:
		return False
	locs = get_locations(iloc, loc_file)
	for loc in locs:
		res = check_loc(loc, peptide, fasta, length_dict)
		if res:
			return res
	return False

# get all matching data for a peptide in a list (desc, seq)
def get_peptide_matches(peptide, fasta, loc_file, ii_dict, length_dict):
	key = peptide[:5]
	if key in ii_dict:
		iloc = int(ii_dict[peptide[:5]])
	else:
		return []
	locs = get_locations(iloc, loc_file)
	out_list = []
	for loc in locs:
		res = check_loc(loc, peptide, fasta, length_dict)
		if res:
			out_list.append(res)
	return out_list


def load_indexes(fasta_name):
	name = '.'.join(fasta_name.split(".")[:-1])
	fasta = open(fasta_name,'rU')
	print "Loading Fasta Index"
	length_dict = cPickle.load(open(name + ".lin"))
	print "Length Index Loaded"
	ii_dict_file = open(name + ".iin")
	loc_dict_file = open(name + ".sin")
	ii_dict = {}
	for line in ii_dict_file:
	    line = line.split(",")
	    if len(line) > 1:
	        ii_dict[line[0]] = line[1].strip()
	print "Secondary Index Loaded"

	return fasta, loc_dict_file, ii_dict, length_dict


'''
Creating
'''

def read_until_EOL(rfile):
	return rfile.readline()
def export_dict(tdict, lfile, iifile):
	for k,v in tdict.items():
		iifile.write("%s,%s\n" % (k, lfile.tell()))
		lfile.write("%s" % k)
		for l in v:
			lfile.write(",%s" % l)
		lfile.write("#")

def create_fasta_index(fasta, word_size):
	file_name = '.'.join(fasta.split(".")[:-1])
	fasta = open(fasta,'rU')
	index_file = open(file_name + ".sin", 'w')
	index_index_file = open(file_name + ".iin", 'w')
	length_index_file = open(file_name + ".lin", 'w')

	seq = ''
	char = 'go'
	cursor_start, cursor_length = 0,0
	index_dict = {}
	length_dict = {}
	while char != '':
		char = fasta.read(1)
		if char == '>' or char == '':
			# throw in last sequence
			index_dict.setdefault(seq,[]).append(cursor_start)
			cursor_length = fasta.tell() - cursor_start -2 #account for new line and character (>,'')
			length_dict[cursor_start] = cursor_length
			read_until_EOL(fasta)
			cursor_start = fasta.tell()
			
			#print cursor_start,cursor_length
			seq = ''
		else:
			if len(seq) < word_size:
				seq += char
			elif char != '\n':
				index_dict.setdefault(seq,[]).append(cursor_start)
				seq = seq[1:] + char
				
	index_dict.pop('',None)
	length_dict.pop(0,None)
	#cPickle.dump(index_dict,index_file, protocol=-1)
	export_dict(index_dict,index_file, index_index_file)
	cPickle.dump(length_dict,length_index_file, protocol=-1)
	index_file.close()
	index_index_file.close()
	length_index_file.close()
	fasta.close()
	return True

def create_tag_index(transcripts_fasta, tag_regex):
	print "Indexing Tags"
	fasta = open(transcripts_fasta, 'rU')
	tag_dict = {}
	char = 'go'
	while char != '':
		char = fasta.read(1)
		if char == '>':
			loc = fasta.tell() - 1 
			test = read_until_EOL(fasta)
			match = tag_regex.search(test)
			if match != None:
				tag = match.groups()[0]
				tag_dict[tag] = loc
	print "Indexing Complete, %i Tags Found" % len(tag_dict)
	fasta.close()
	return tag_dict




