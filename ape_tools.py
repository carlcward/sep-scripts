# convert a codon index to an dna index thats starts from 1
def index_frame_to_loc(index,frame):
    return (index*3 + frame + 1)

# and convert a ape location back to an index, frame
def loc_to_index_frame(loc):
    frame = (loc-1) % 3
    index = ((loc-1) - frame) / 3
    return index,frame

# calculate the location in the protein based on the ape locations for the protein
def calculate_location_in_protein(start, stop, cds_start, cds_stop):
    location = ''
    if cds_start == 'non-coding':
        location = 'non-coding'
        cds_start = 0
        cds_stop = 0
    elif start < loc_to_index_frame(int(cds_start))[0]:
        
        location = "5' UTR"
    elif stop > loc_to_index_frame(int(cds_stop))[0]:
        location = "3' UTR"
    else:
        
        location = 'CDS'
    return location, cds_start, cds_stop

# create the string to codes for the ape file
def create_ape_map(dna, frame, p_st, p_sp, _start, _stop, p_length, cds_start, cds_stop):
	# ape file indexing is INCLUSIVE 
	# so 1..3 highlights 1, 2, and 3 
    length = len(dna)
    p_start = index_frame_to_loc(p_st, frame)
    p_stop = index_frame_to_loc(p_sp, frame)-1
    start = index_frame_to_loc(_start, frame)
    stop = index_frame_to_loc(_stop, frame)

    ape_string = '''LOCUS       New_DNA                 %i bp ds-DNA     linear       17-JUN-2013
DEFINITION  .
ACCESSION   
VERSION     
SOURCE      .
  ORGANISM  .
COMMENT     auto-generated by python script from given protein and sequence, calculated protein length = %i
COMMENT		adsfasdf
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    %s..%s
                     /label=CDS
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    %i..%i
                     /label=New Feature
                     /ApEinfo_fwdcolor=#ff8000
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    %i..%i
                     /label=New Feature(1)
                     /ApEinfo_label=New Feature
                     /ApEinfo_fwdcolor=#ff0000
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    %i..%i
                     /label=New Feature(2)
                     /ApEinfo_label=New Feature
                     /ApEinfo_fwdcolor=#ff0000
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
      %s 
//
''' % (length, p_length, cds_start, cds_stop, p_start, p_stop, start, start+2, stop, stop+2, dna)
    return ape_string