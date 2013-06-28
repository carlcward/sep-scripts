#!/usr/bin/python
import argparse, StringIO
import xml.etree.ElementTree as etree
from urllib import urlopen, urlencode
from time import sleep
import re
#set up the parsers
#parser = argparse.ArgumentParser('Online Blast Search')
#parser.add_argument("--fasta", "-f",dest="ref", type=str)
#parser.add_argument("--pep", "-p", dest="pep", type=str)
#parser.add_argument("--out", "-o", dest="out", type=str)
#args = parser.parse_args()

'''
ref_file = open(args.ref, "r")
pep_file = open(args.pep, "r")

out_file = open(args.out, "w")
'''
delchars = ''.join(c for c in map(chr, range(256)) if not c.isalpha())
def send_blast_request(queries):
	word_size = 2
	expect_value = 5000000
	matrix_name = "PAM30"
	database = "refseq_rna"
	#query_file = pep_file.read()

	querystring = ""
	for q in queries:
		querystring += (">%s\n%s\n" % (q,q))

	data = {
		'QUERY' : querystring,
		'DATABASE' : database,
		'WORD_SIZE' : word_size,
		'EXPECT' : expect_value,
		'MATRIX' : matrix_name,
		'HITLIST_SIZE' : 100,
		'PROGRAM' : 'tblastn',
		'NCBI_GI' : 'on',
		#'PAGE' : 'Proteins',
		'ENTREZ_QUERY' : 'Homo sapiens[Organism]',
		'CMD' : 'Put'
	}
	url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?" + urlencode(data)
	print "Sending Blast Query for %i Peptides\n" % len(queries)
	return urlopen(url)
def get_blast_results(put_request):
	rid = None
	rtoe = None
	for line in put_request:
		if line[:7] == "    RID":
			rid = line.split("=")[1].strip()
		elif line[:8] == "    RTOE":
			rtoe = int(line.split("=")[1].strip())

	get_request = urlopen(("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?" +
	        "CMD=Get&RID=%s&FORMAT_OBJECT=Alignment"+
	        "&FORMAT_TYPE=XML&DESCRIPTIONS=100&ALIGNMENTS=200" +
	        "&ALIGNMENT_TYPE=Pairwise") % (rid))

	print "Waiting for Blast results\nEstimate wait time %i seconds\n" % (rtoe*2)
	sleep(rtoe/2)
	done = False
	while not done:
		data = get_request.read()
		#print data
		if re.search("xml version=\"1.0\"", data):
				done = True

		else:

			sleep(rtoe/3)
			get_request = urlopen(("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?" +
	        "CMD=Get&RID=%s&FORMAT_OBJECT=Alignment"+
	        "&FORMAT_TYPE=XML&DESCRIPTIONS=100&ALIGNMENTS=200" +
	        "&ALIGNMENT_VIEW=Pairwise") % (rid))
	print "Blast results received"
	return data
def parse_blast_results(xml):
	xml_file = open("blast_xml.txt", 'w')
	xml_file.write(xml)
	xml = StringIO.StringIO(xml)

	tree = etree.parse(xml)
	root = tree.getroot()
	queries = root.find('BlastOutput_iterations')
	query_dict = {}

	print "Processing Blast Results"

	for query in queries:
		# iterate through the peptides with any matches
		peptide = query.find('Iteration_query-def').text
		length = int(query.find('Iteration_query-len').text)
		hit_list = []
		for hit in query.find('Iteration_hits'):
			# get the data for each hit
			hit_id = hit.find('Hit_id').text
			hit_val = int(hit_id.split("|")[1])
			hit_name = hit.find('Hit_def').text.split('|')[0]
			# filter out ONLY the perfect matches (indentities = peptide length)
			hit_hsps = filter(lambda h: h.find('Hsp_identity').text == str(length), hit.find('Hit_hsps'))

			if len(hit_hsps) > 0:
				#print "HIT: %s\t%s" % (peptide,hit_name)

				# fetch the data for the gene from genbank
				ncbi_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=%i&db=nuccore&dopt=genbank&retmode=text" % (hit_val)
				record = urlopen(ncbi_url).read()
				# look only through the features
				record = record.split("FEATURES")[1] 
				record  = StringIO.StringIO(record)
				cds_range = ''
				seq = ''

				# get the coding sequence location and mRNA sequence, return it with the names
				for line in record:
					cds_search = re.search("\WCDS\W(.*)",line)
					if cds_search:
						cds_range = cds_search.group(1).strip()
					if 'ORIGIN' in line:
						while not line.startswith("//"):
							line = record.next()
							seq += line.translate(None,delchars).upper()
				# append name, coding range, and sequence to the hitlist for each query
				# DO NOT APPEND IF NO CDS 
				if cds_range.strip() != '':
					hit_list.append((hit_name, hit_id, cds_range, seq))
				else:
					hit_list.append((hit_name, hit_id, 'non-coding..non-coding', seq))
		if len(hit_list) > 0:
			query_dict[peptide] = hit_list
	return query_dict



if __name__ == '__main__':
	p = send_blast_request(["AGSGKNDHAEKV","AMRYVASYLLAALGGNSSPSAKD","EVISELNGKN"])
	xml = get_blast_results(p)
	parse_blast_results(xml)
