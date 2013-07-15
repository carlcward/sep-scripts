import re, sys
from collections import defaultdict
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
			code_string = code_dict[code][0] + " - " + code_dict[code][1]
			coord_string = "chr%s:%i-%i" % (chrom, coord_start, coord_end)
			coordinate_map[coord_key].append((coord_start,coord_end, [coord_string] +[line[0]] +[code_string] + line[7:16]))
	return coordinate_map

def read_sep_file(file_name):
	coordinate_map = defaultdict(list)
	sep = open(file_name, 'rU')
	for line in sep:
		line = line.split("\t")
		#print line
		tcoords = line[7].split(";")
		tcoords = filter(lambda c: c.strip() != '', tcoords)
		for coord_string in tcoords:
			regex = re.match("chr([0-9]{1,2}):([0-9]*)\-([0-9]*)",coord_string)
			if regex != None:
				# offset of sep start from transcript coords, sep length in bases
				# transcript and transcript coords do not match up however... because of various isofroms
				# so for now just looking at entire trasncrupt coordinates
				#sep_start_offset = int(line[7]) - 1
				#sep_length_bases = int(line[5])*3
				regex = regex.groups()
				coord_key = regex[0]
				# sep coord start = transcript coord + sep offest
				coord_start = int(regex[1]) #+ sep_start_offset
				coord_end = int(regex[2]) #coord_start + sep_length_bases
				#print coord_start,coord_end
				coordinate_map[coord_key].append((coord_start,coord_end,line[:6]))
	return coordinate_map

if __name__ == '__main__':
	# collect the matches... where a known mutation lies within a TRANSCRIPT for one of the SEPs 
	# (not necessarily within the sep itsself). We need a way to start the exact location info for the SEPs 
	hits = []
	code_dict = create_tissue_source_dict(TISSUE_SOURCE_FILE_PATH)
	sep_coord = read_sep_file(sys.argv[1])
	maf_coord = read_maf(sys.argv[2], code_dict)
	# print the lengths of each dict (all values) (fancy flatten syntax)
	print (len([item for sub in sep_coord.values() for item in sub]),len([item for sub in maf_coord.values() for item in sub]))
	# for chromosome with seps 
	for chrom,data in sep_coord.items():
		# for each sep on a transcript
		for sep_start,sep_stop,sep_data in data:
			# for each mutation on the same chromosome
			for mut_start,mut_stop,mut_data in maf_coord[chrom]:
				# see if the mutation matches with the sep
				if mut_start > sep_start and mut_start < sep_stop:
					hits.append((sep_data,mut_data))
	# collected results
	print len(hits)
	for sep,mut in hits:
		print "%s\n%s\n" % ('\t'.join(sep),'\t'.join(mut))

