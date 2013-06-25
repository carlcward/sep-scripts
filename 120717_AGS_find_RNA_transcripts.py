from fasta_reader import Data_Reader
import sys
alphabet=['A',"B","C","D",'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
plist=[]

# syntax: peptide_file translated_database transscripts coord_database 

def PCnumber(match):
    #end=match.find('_',4)
    PCtag=match.split('-')[0]
    print "Tag: " + PCtag
    return PCtag

#write=open("PCtags.txt","a")

## searches database for PCtags
def Searchtranslateddatabase(peptide):
    out_list = []
    infile=open(sys.argv[2],'rU')
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
    

## searches translated database for peptides in the supplied file, and collects answers in a list
PCtags = []
query=Data_Reader(sys.argv[1])
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]
pqeueries = []
for line in query:
    pqeueries.append(line)
for line in f7(pqeueries):
    print "INPUT: " + line
    PCtags += Searchtranslateddatabase(line)

##closes all opened files
query.close()
##This section matches PC# vs TCON#

##searches RNA database with PCnumber, returns list of matching transscripts
def SearchRNA(input,previous_tag):
    trans_list = []
    RNAfile=open(sys.argv[3],'r')
    print "Searching!"
    print input
    for transcripts in RNAfile:
        #if transcripts.find("_"+input[1:])!=-1:
        if 'TCONS_' + input in transcripts:
            print "Found " + input 
            transcripts=RNAfile.next()
            assembly=""
            while transcripts.find(">")==-1:
                assembly=assembly+transcripts[0:-1]
                transcripts=RNAfile.next()
            trans_list.append(assembly)
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

#transcriptsfile=open("tideswithtranscripts.txt","w")
data_file = open("peptides_with_data.txt","w")


# put all the coords in a dictionary with key by PCNumber
coord_file = open(sys.argv[4])

print "Fetching Coordinates..."
coord_dict = get_coordinates_from_DTASelect(coord_file)
'''
for line in coord_file:
    data = line.split("\t")
    # look for the transcript lines
    if 'transcript' in data[2]:
        #range=chr9:127023982-127024178 strand=+
        # extract PCnum and coord
        pcnum = data[8].split('"')[3]
        coord_string = "%s:%s-%s strand=%s" % (data[0],data[3],data[4],data[6])
        coord_dict[pcnum] = coord_string
'''


# go through all the stored peptide/pctags and write the data out
# matches PCtags with the PCtags on the RNA transcripts
for peptide,tag in PCtags:
    transcripts = SearchRNA(tag,peptide)  
    for transcript in transcripts:
        data_file.write("%s\t%s\t%s\t%s\n" % (coord_dict.get(tag,'COULD NOT FIND'),tag,peptide,transcript))


data_file.close()

