from fasta_reader import Data_Reader
from fa_index import get_locations, get_peptide_matches, load_indexes, create_fasta_index, create_tag_index
import sys, cPickle, re
plist=[]

# syntax: peptide_file translated_database transscripts coord_database 

# get the right regex for seraching through the transcript file and creating the tag index
def get_matching_transcript_file_regex(test):
    possible_tag_regexs = ( x for x in [
    "([a-z0-9]*\.[0-9])",
    "\|((?:AC|AP|NC|NG|NM|NP|NR|NT|NW|NZ|XM|XP|XR|YP|ZP)_.*?)\|",
    "TCONS_([0-9]*)"]
    )
    
    regex = None
    try:
        while not regex:
            regex = re.search(possible_tag_regexs.next(), test)
        matching_regex = regex.re.pattern
        
        return matching_regex
    except:
        print "None!"
        return None

# return the PC number from the given info line from the translated fasta file
# automatically detectcs the database tag format
def PCnumber(match):
    possible_tag_regexs = ( x for x in [
    "^>([a-z0-9]*\.[0-9])",
    "^>((?:AC|AP|NC|NG|NM|NP|NR|NT|NW|NZ|XM|XP|XR|YP|ZP)_[0-9]*\.[0-9])[a-z]",
    "^>PC_([0-9]*)_"]
    )
    regex = None
    print match
    try:
        while not regex:
            regex = re.search(possible_tag_regexs.next(), match)
        PCtag = regex.groups()[0]
        #PCtag = match.split('_')[1]
        print ("Tag: " + PCtag)
        return PCtag
    except:
        print "Unrecognized Tag Format"
        return None

#write=open("PCtags.txt","a")

## searches database for PCtags
'''
def Searchtranslateddatabase(peptide):
    out_list = []
    infile=cPickle.load(sys.argv[2])

    previous_line=""
    up_one_line=""
    
    for lines in infile:
        if (up_one_line[:-1]+lines).find(peptide)!=-1:
                if (peptide[:-1]+previous_line) not in plist:
                    #write.write(peptide+"\n")
                    #write.write(">"+str(PCnumber(previous_line)))
                    #write.write('\n')
                    out_list.append((peptide,PCnumber(previous_line)))
                    #print up_one_line+lines
                    plist.append(peptide[:-1]+previous_line)
        up_one_line=lines
        if lines.find(">")!=-1:    
            previous_line=lines
    return out_list        
    infile.close()
'''
def Searchtranslateddatabase(peptide, ii_dict, length_dict, loc_file, ffile):
    # get the (desc,seq) of the matching sequences and create the necessary list of (pep, PCtag)
    return map(lambda a: (peptide, PCnumber(a[0])),get_peptide_matches(peptide, ffile, loc_file, ii_dict, length_dict))
    '''
    iloc = int(ii_dict[peptide[:5]])
   
    locs = get_locations(iloc, loc_file)
    
    out_list = []
   
    
    for loc in locs:
        
        
        seq, desc = check_loc(loc, peptide, ffile, length_dict)
        out_list.append((peptide,PCnumber(desc[::-1])))

    return out_list
    '''



    


##This section matches PC# vs TCON#

##searches RNA database with PCnumber, returns list of matching transscripts
def SearchRNA(input,previous_tag, tag_dict):
    trans_list = []
    RNAfile=open(sys.argv[3],'r')
    print ("Searching!")
    #print (input)
    
    try:
        RNAfile.seek(tag_dict[input])
        transcripts = RNAfile.readline()
        if re.search("[^0-9]"+input+"[^0-9]",transcripts):
            print ("Found " + input )
            transcripts=RNAfile.next()
            assembly=""
            while transcripts.find(">")==-1:
                assembly=assembly+transcripts[0:-1]
                transcripts=RNAfile.next()
            trans_list.append(assembly)
    except KeyError:
        print 'Could Not Find ' + input
        pass
    return trans_list
                
            
def get_coordinates_from_DTASelect(coord_file):
    out_dict = {}
    for line in coord_file:
        data = line.split("\t")
        #print len(data),data[0]
        if data[0][:2] == 'PC':
            PC = data[0].split('_')[1]
            for datum in data:
                if datum.startswith('range'):
                    coord = datum.split('RNAsource')[0].strip()
                    out_dict[PC] = coord
    return out_dict

def get_coordinates_from_transcripts_gtf(gtf_file):
    out_dict = {}
    for line in gtf_file:
        data = line.split("\t")
        # look for the transcript lines
        if 'transcript' in data[2]:
            #range=chr9:127023982-127024178 strand=+
            # extract PCnum and coord
            pcnum = data[8].split('"')[3]
            coord_string = "%s:%s-%s strand=%s" % (data[0],data[3],data[4],data[6])
            out_dict.setdefault(pcnum, []).append(coord_string)
        if 'exon' in data[2]:
            pcnum = data[8].split('"')[3]
            exon_start = int(data[3])
            exon_end = int(data[4])
            out_dict.setdefault(pcnum, []).append((exon_start,exon_end))
    return out_dict

# go through all the stored peptide/pctags and write the data out
# matches PCtags with the PCtags on the RNA transcripts
def write_data_file(data_file, PCtags, coord_dict, tag_dict):
    for peptide,tag in PCtags:
        transcripts = SearchRNA(tag,peptide, tag_dict)  
        for transcript in transcripts:
            out_string = "%s\t%s\t%s\t%s\t \n" % (coord_dict.get(tag,['No Coord'])[0],tag,peptide,transcript.strip())
            data_file.write(out_string)

# write the exact same outfile, except with the exon data gleaned fromt he GTF file
# format = 
def write_data_file_with_exons(data_file, PCtags, coord_dict, tag_dict):
    for peptide,tag in PCtags:
        transcripts = SearchRNA(tag,peptide,tag_dict)  
        for transcript in transcripts:
            data_file.write("%s\t%s\t%s\t%s\t" % (coord_dict.get(tag,['No Coord'])[0],tag,peptide,transcript))
            for es,ee in coord_dict.get(tag,[])[1:]:
                data_file.write("%i,%i;" % (es,ee))
            data_file.write(" \n")



## searches translated database for peptides in the supplied file, and collects answers in a list
PCtags = []

try:
    trans_file,loc_dict_file,ii_dict,length_dict = load_indexes(sys.argv[2])
except IOError:
    print "Index Not Found\nCreating Index - This requires ~10min and significant amounts of RAM\n"
    create_fasta_index(sys.argv[2],5)
    print "Indexing Complete"
    trans_file,loc_dict_file,ii_dict,length_dict = load_indexes(sys.argv[2])

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


pqeueries = []
query=Data_Reader(sys.argv[1])
for line in query:
    pqeueries.append(line)
count = 0

# remove duplicates but stay in order...
for line in f7(pqeueries):
    count += 1
    if count % 100 == 0:
        print "%i Peptides Done!" % count
    print ("INPUT: " + line)
    #fetch all the tags
    PCtags += Searchtranslateddatabase(line,ii_dict,length_dict,loc_dict_file,trans_file)

transcript_file = open(sys.argv[3])
# get the regex that will extract the tags from the transcript file
tag_regex = get_matching_transcript_file_regex(transcript_file.readline().strip())
transcript_file.close()
tag_dict = create_tag_index(sys.argv[3], re.compile(tag_regex))
trans_file.close()
##closes all opened files
query.close()

data_file = open("peptides_with_data.txt","w")


# put all the coords in a dictionary with key by PCNumber
coord_file = open(sys.argv[4])

print ("Fetching Coordinates...")

#DTA SELECT VS GTF FILE
#coord_dict = get_coordinates_from_DTASelect(coord_file)
coord_dict = get_coordinates_from_transcripts_gtf(coord_file)

#write_data_file_with_exons(data_file, PCtags, coord_dict, tag_dict)
write_data_file_with_exons(data_file,PCtags,coord_dict, tag_dict)

data_file.close()
trans_file.close()
loc_dict_file.close()

