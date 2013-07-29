import re, sys
from collections import defaultdict, OrderedDict
'''
Reads in mutation files from TCGA data and checks to see if any of the seps has a mutation in item
Required addition of additional output from sep_processor - where in the transcript each sep start (base offset)
so the coordinates of the mutation and sep could be lined up.
'''
TISSUE_SOURCE_FILE_PATH = '/Users/carlcward/Desktop/Carl/MAFs/tissueSourceSite.txt'
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def create_tissue_source_dict(tss_file_name):
	tss_file = open(tss_file_name,'rU')
	out_dict = {}
	for line in tss_file:
		line = line.split("\t")
		code = line[0].strip()
		name = line[2].strip()
		short = line[3].strip()
		out_dict[code] = (short,name)
	tss_file.close()
	return out_dict


def read_maf(file_name, code_dict):
	coordinate_map = defaultdict(list)
	maf = open(file_name, 'rU')
	for line in maf:
		line = line.split("\t")
		
		if len(line) > 5 and RepresentsInt(line[5]):
			# chromosome, mutation start coord, mutation stop coord
			
			chrom = line[4]
			coord_start = int(line[5])
			coord_end = int(line[6])
			coord_key = "%s" % (chrom)
			code = line[16].split("-")[1]
			code_string = code_dict[code][1]
			coord_string = "chr%s:%i-%i" % (chrom, coord_start, coord_end)
			coordinate_map[coord_key].append((coord_start,coord_end, [coord_string] +[line[0]] +[code_string] + line[7:16]))
	maf.close()
	return coordinate_map

# get the index of sep transcipts coords in the sep_processor results file (dynamic)
def read_header(sep_file):
	header = sep_file.readline().split("\t")
	index = 0
	for i in range(len(header)):

		if  'Sep Start In Transcript' in header[i].strip() or 'Sep Transcript Coords' in header[i].strip():
			return i


def read_sep_file(file_name):
	coordinate_map = OrderedDict()
	sep = open(file_name, 'rU')
	sep_coords_i = read_header(sep)
	for line in sep:
		line = line.split("\t")
		peptide_data = line[:9]
		peptide = line[1]
		peptide_key = tuple(peptide_data[8])
		if peptide_key not in coordinate_map:
			coordinate_map[peptide_key] = [peptide_data]
			tcoords = line[sep_coords_i].split(";")
			tcoords = filter(lambda c: c.strip() != '', tcoords)
			for coord_string in tcoords:
				# transcript coords from the Sep Start Column in out from sep_processor.py
				regex = re.match("chr([0-9XY]{1,2}):([0-9]*)\-([0-9]*)",coord_string)
				if regex != None:
					
					
					regex = regex.groups()
					coord_key = regex[0]
					# sep coord start = transcript coord + sep offest
					coord_start = int(regex[1]) #+ sep_start_offset
					coord_end = int(regex[2]) #coord_start + sep_length_bases
					# the exact same data already exist in the dictions, do not add it
					# other wise add in the form (start coord, end coord, peptide data)

					coordinate_map[peptide_key].append((coord_key,coord_start,coord_end))
	sep.close()
	return coordinate_map

if __name__ == '__main__':
	# collect the matches... where a known mutation lies within a TRANSCRIPT for one of the SEPs 
	# (not necessarily within the sep itsself). We need a way to start the exact location info for the SEPs 
	hits = OrderedDict()
	code_dict = create_tissue_source_dict(TISSUE_SOURCE_FILE_PATH)
	sep_coord = read_sep_file(sys.argv[1])
	maf_coord = read_maf(sys.argv[2], code_dict)
	out_file = open(sys.argv[3], 'w')
	# print the lengths of each dict (all values) (fancy flatten syntax)
	print (len([item for sub in sep_coord.values() for item in sub]),len([item for sub in maf_coord.values() for item in sub]))
	# for unique sep with seps 
	for sep,data in sep_coord.items():
		sep_data = data[0]
		loc = 1
		# for each sep on a transcript
		
		for chrom,sep_start,sep_stop in data[1:]:
			#print sep_stop,sep_start, sep
			
			# for each mutation on the same chromosome
			for mut_start,mut_stop,mut_data in maf_coord[chrom]:
				# see if the mutation matches with the sep but dont include the same mutation for the same set of sep data
				# tuples are hashable to create a dict with them of sep data for duplicate exclusion
				if mut_start > sep_start and mut_start < sep_stop:
					#print loc, mut_start-sep_start
					mut_id = "%s,%s" % (loc + mut_start-sep_start, mut_stop-mut_start)
					
					if mut_id not in [mut[0] for mut in hits.get(tuple(sep_data),[])]:
						hits.setdefault(tuple(sep_data),[]).append([mut_id]+mut_data)
			loc += sep_stop - sep_start
	# collected results
	print len(hits)
	for sep,v in hits.items():
		outfile.write('\n')
		for mut in v:
			outfile.write("%s\n%s\n" % ('\t'.join(('',)+sep),'\t'.join(mut)))
	outfile.close()


