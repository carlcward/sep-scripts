import os
'''
Reads FASTA files return description and sequence, or just sequence
'''
class Fasta_Reader:
	def __init__(self, file_path):
		self.fasta = open(file_path, "rU", os.O_NONBLOCK)
		self.prev_desc = ''
	

	# the class can be used as an iterator
	def __iter__(self):
		data = self.next_data()
		while data != None:
			yield data
			data = self.next_data()

	
	
	def next_data(self):
		seq = ''
		line = self.fasta.readline()
		while True:
			if line == '' or line[0] == '>' :
				out = (self.prev_desc[:],seq)
				self.prev_desc = line.strip()
				if seq != '':
					return out
				if line == '':
					return None
			else:
				seq += line.strip()
			line = self.fasta.readline()

	def next_sequence(self):
		return self.next_data()[1]
	def close(self):
		self.fasta.close()
'''
if __name__  == '__main__':
	FF = Fasta_Reader("../myFasta.fa")
	print FF.next_sequence()
	print FF.next_data()
	print "####################"
	for data in FF:
		print data
'''

'''
Can read both line by line AND FASTA files
Will skip over description lines, and work with windows delimited lines
'''

class Data_Reader:
	
	delchars = ''.join(c for c in map(chr, range(256)) if not c.isalpha())
	def __init__(self, file_path):
		self.pfile = open(file_path, 'rU', os.O_NONBLOCK)
		

	def __iter__(self):
		seq = self.next_sequence()
		while seq != None:
			yield seq
			seq = self.next_sequence()

	def next_sequence(self):
		line = self.pfile.readline()
		if line == '':
			return None
		elif line.strip() == '' or line.startswith('>'):
			return self.next_sequence()
		else:
			return line.strip().translate(None,self.delchars)
	def close(self):
		self.pfile.close()

		

if __name__ == '__main__':
	dd = Data_Reader('../database_test/nblast_in.txt')
	for line in dd:
		print line.strip()