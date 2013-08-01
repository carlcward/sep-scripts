import sys, subprocess, os
import twill.commands as tc
import Image
import pyPdf
from fasta_reader import Data_Reader
#from screenshot import Screenshot

# The base site url for the ip2
SITE_URL = "http://ip2-saglab.fas.harvard.edu"

# File to load peptides from
pep_file = Data_Reader(sys.argv[1])

# get the url of the 'spreadsheet' from copy paste -> stdin
print "Copy and pasta the URL of the protein page from your web browser:"
data_file = sys.stdin.readline().strip()
data_file += "&display=peptide"

# the path to save the spectra
out_path = sys.argv[2]
out_path = out_path.rstrip('/') + "/"
# make the outpath if it doesnt exist
if os.path.isdir(out_path):
	print "Out directory exists"
else:
	print "Out directory does not exist, Created!"
	subprocess.call(['mkdir', out_path])

# collect the peptides (can take FASTA or line by line peptides because of set)
peps = set()

# quick method for stripping and deleting ALL but letters
delchars = ''.join(c for c in map(chr, range(256)) if not c.isalpha())
for line in pep_file:
	peps.add(line)
# clean up peps
peps = filter(lambda p: p.strip() != '' , peps)
num_starting_peps = len(peps)

# open our spreadsheet
tc.go(data_file)
print "Navigation Successful"
# login
print "Logging In..."
tc.fv("1", "j_username", "jiao")
tc.fv("1", "j_password", "jiao321")
tc.submit('4')
links = []
print "Login Successful, Fetching peptides..."
for l in  tc.showlinks():
	links.append(l.url)


links = filter(lambda l: 'Lorikeet' in l,links)

# method to grab the sequence from the URL
def extract_seq(l):
	seq = l.split('sequence')[1].split('&')[0][3:-2]
	seq = seq.translate(None, delchars)
	return seq

# filter method that keeps track of which peps have been used
def in_pep(l):
	
	seq = extract_seq(l)
	
	if seq in peps:
		peps.remove(seq)
		return True
	else:
		return False


links = filter(in_pep, links)


'''
SCREENGRABBING
'''
# java script to log in again

count = 0;
filenames = []
for link in links:
	js = '''
	document.getElementById('j_username').value = 'jiao';
	document.getElementById('j_password').value = 'jiao321';
	document.getElementById('loginForm').submit();
	''' 
	count += 1
	# call the screen grabber
	subprocess.call([
		'webkit2png',
		'--delay=.75',
		'-F',
		'--js=' + js,
		SITE_URL + link.replace("chargeState=2","chargeState=3"),
		'--dir=' + out_path,
		'--filename=' + extract_seq(link)
		],stdout=open(os.devnull, 'wb'))
	filenames.append(extract_seq(link))
	print "\n%i of %i Images Collected\n" % (count, len(links))


print "Image Post-processing..."
# crop the images and delete the old ones (370 from top 85 from bottom)
for f in filenames:
	img = Image.open(out_path+f+'-full.png')
	subprocess.call(['rm',out_path+f+'-full.png'])
	width,height = img.size
	box = (22,370,width,height-85)
	final_img = img.crop(box)
	out_file = open(out_path+f + '.png','w')
	final_img.save(out_file,'PNG')
	out_file.close()

# create the PDFs

#out_pdf = open(outpath + "spectra.pdf", "w")
#for f in filenames:


print "Input File had %i Peptides" % num_starting_peps
print "%i Matching Links Found" % len(links)
if len(peps) > 0:
	print "Complete... Failed to Find %i Peptides" % len(peps)
	print peps
else:
	print "Complete... All Peptides Found"
