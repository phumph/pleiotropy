import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re 
from Bio import SeqIO

#define working directory from which this program is being run
wkdir = '/Volumes/Promise\ Pegasus/Lucas/whole_genome_seq/CLM_clones/'
#wkdir = '~/Desktop/SherlockLab/revision_work/sequence_analysis/'

#define directory where all the fastq data files are and file prefix and suffix
fastq_Dir = 'fastq_files/'
#fastqPrefix = '141126_PINKERTON_0343_BC4J1PACXX_L7_'
r1FileSuffix = '_R1_001.fastq'
r2FileSuffix = '_R2_001.fastq'

get_BC_in_fastq = False
most_freq_BC = True

h_freq_file = open('most_freq_bc.txt', 'w')

#all_index_tags_dic  = {'P1-D10':'ACTGATCG-GCGTAAGA_S46', 'P1-D1':'TAGCGAGT-GCGTAAGA_S37'}
all_index_tags_dic  = {}
index_file = open ('index_key_NextSeq_DY_BFA_rev.txt', 'r')
sample_count = -1
for line in index_file:
	(strain, I7, seq7, I5, seq5) = line.rstrip('\n').split('\t') # split the line
	sample_count += 1
	if strain.startswith('P') :
		all_index_tags_dic[strain] = seq7 + '-' + seq5 + '_S' + str(sample_count)
index_file.close()


for strain in all_index_tags_dic:
	index_tag = strain + '-' + all_index_tags_dic[strain]
	print (index_tag)
	#TAG SPECIFIC PARAMETERS
	#define the index tag data to be analyzed
	#index_tag = 'TAAGGCGA-TAGATCGC'

	#define files names of this particular index tag
	fastqR1 = fastq_Dir + index_tag + r1FileSuffix
	fastqR2 = fastq_Dir + index_tag + r2FileSuffix

	forwardFQIterator = SeqIO.parse(fastqR1,"fastq")
	reverseFQIterator = SeqIO.parse(fastqR2,"fastq")

	#make a folder for data output for this clone
	outputDir = index_tag + "_data/"
	
	BC_both_dir = outputDir + 'Extracted_BC.txt'
	count_both_dir = outputDir + 'Count_both_bc.txt'
	count_bc1_dir = outputDir + 'Count_bc1_bc.txt'
	count_bc2_dir = outputDir + 'Count_bc2_bc.txt'
	count_bc_dir = outputDir + 'Count_all_bc.txt'
	
	

	if get_BC_in_fastq == True:

	
		#pattern_BC1 = re.compile('(.GTTAT|A.TTAT|AG.TAT|AGT.AT|AGTT.T|AGTTA.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.GTACC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)')
		#pattern_BC2 = re.compile('(.GTACC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)\D{4,7}?AA\D{4,7}?AA\D{4,7}?TT\D{4,7}?(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC.)')

			
		#pattern_rev_BC1 = re.compile('(.GTACC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC.)')
		#pattern_rev_BC2 = re.compile('(.GTTAT|A.TTAT|AG.TAT|AGT.AT|AGTT.T|AGTTA.)\D{4,7}?AA\D{4,7}?TT\D{4,7}?TT\D{4,7}?(.GTACC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)')

		pattern_BC1 = re.compile('(ATACGAAGTTAT)\D{5}?AA\D{5}?AA\D{5}?TT\D{5}?(GGTACCGATATC)')
		pattern_BC2 = re.compile('(GATATCGGTACC)\D{5}?AA\D{5}?AA\D{5}?TT\D{5}?(ATAACTTCGTAT)')

			
		pattern_rev_BC1 = re.compile('(GATATCGGTACC)\D{5}?AA\D{5}?TT\D{5}?TT\D{5}?(ATAACTTCGTAT)')
		pattern_rev_BC2 = re.compile('(ATACGAAGTTAT)\D{5}?AA\D{5}?TT\D{5}?TT\D{5}?(GGTACCGATATC)')


	
		bc1_dic = {}
		bc2_dic = {}
		both_dic = {}
	


		for fwd_record in forwardFQIterator :
			forwardSeq = str(fwd_record.seq)
			seq_id = str(fwd_record.id)
			BC1_match = re.search(pattern_BC1, forwardSeq)
			BC2_match = re.search(pattern_BC2, forwardSeq)
	
			if BC1_match and BC2_match :
				BC1_seq = BC1_match.group() # take the sequence that match the pattern
				BC2_seq = BC2_match.group()
				BC1 = BC1_seq[12:len(BC1_seq)-12]
				BC2 = BC2_seq[12:len(BC2_seq)-12]
				both_bc = BC2 + '_' + BC1
				both_dic[seq_id] = both_bc

			
			elif BC2_match:
				BC2_seq = BC2_match.group()
				BC2 = BC2_seq[12:len(BC2_seq)-12]
				bc2_dic[seq_id] = BC2

		
			elif BC1_match:
				BC1_seq = BC1_match.group()
				BC1 = BC1_seq[12:len(BC1_seq)-12]
				bc1_dic[seq_id] = BC1


			
		for rev_record in reverseFQIterator :
			reverseSeq = str(rev_record.seq)
			seq_id_rev = str(rev_record.id)
			BC1_rev_match = re.search(pattern_rev_BC1, reverseSeq)
			BC2_rev_match = re.search(pattern_rev_BC2, reverseSeq)
	
			if BC1_rev_match and BC2_rev_match :
				BC1_rev_seq = BC1_rev_match.group() # take the sequence that match the pattern
				BC2_rev_seq = BC2_rev_match.group()
				BC1_rev = BC1_rev_seq[12:len(BC1_rev_seq)-12]
				BC2_rev = BC2_rev_seq[12:len(BC2_rev_seq)-12]
				mybc1 = Seq(BC1_rev, generic_dna)
				mybc2 = Seq(BC2_rev, generic_dna)
				bc1 = str(mybc1.reverse_complement())
				bc2 = str(mybc2.reverse_complement())
				both_rev_bc = bc2 + '_' + bc1
				both_dic[seq_id_rev + '*' ] = both_rev_bc
			
			elif BC2_rev_match:
				BC2_rev_seq = BC2_rev_match.group()
				BC2_rev = BC2_rev_seq[12:len(BC2_rev_seq)-12]
				mybc2 = Seq(BC2_rev, generic_dna)
				bc2 = str(mybc2.reverse_complement())
				bc2_dic[seq_id_rev + '*' ] = bc2

		
			elif BC1_rev_match :
				BC1_rev_seq = BC1_rev_match.group()
				BC1_rev = BC1_rev_seq[12:len(BC1_rev_seq)-12]
				mybc1 = Seq(BC1_rev, generic_dna)
				bc1 = str(mybc1.reverse_complement())
				bc1_dic[seq_id_rev + '*' ] = bc1
			
	
		count_bc1_dic = {}
		count_bc2_dic = {}
		count_both_dic = {}

		BC_both_file = open (BC_both_dir, 'w')
	
		for seq_id_both in both_dic:
			both_bc = both_dic[seq_id_both]
			bc2, bc1 = both_bc.split('_')
			w_both = seq_id_both + '\t' + bc2 + '\t' + bc1 + '\n'	
			BC_both_file.write(w_both)
		
			if both_dic[seq_id_both] not in count_both_dic:
				count_both_dic[both_bc] = 1
			else:
				count_both_dic[both_bc] += 1
		
		
			if bc1 not in count_bc1_dic:
				count_bc1_dic[bc1] = 1
			else:
				count_bc1_dic[bc1] += 1
			
			
			if bc2 not in count_bc2_dic:
				count_bc2_dic[bc2] = 1
			else:
				count_bc2_dic[bc2] += 1
	
		for seq_id_bc1 in bc1_dic:
			wbc1 = seq_id_bc1 + '\t None' + '\t' + bc1_dic[seq_id_bc1] + '\n'
			BC_both_file.write(wbc1)
		
			if bc1_dic[seq_id_bc1] not in count_bc1_dic:
				count_bc1_dic[bc1_dic[seq_id_bc1]] = 1
			else:
				count_bc1_dic[bc1_dic[seq_id_bc1]] += 1
		
		
		for seq_id_bc2 in bc2_dic:
			wbc2 = seq_id_bc2 + '\t' + bc2_dic[seq_id_bc2] + '\t None\n'		
			BC_both_file.write(wbc2)
		
			if bc2_dic[seq_id_bc2] not in count_bc2_dic:
				count_bc2_dic[bc2_dic[seq_id_bc2]] = 1
			else:
				count_bc2_dic[bc2_dic[seq_id_bc2]] += 1
		
		
		BC_both_file.close()
	
		count_both_file = open(count_both_dir, 'w')
		for key in count_both_dic:
			wboth = str(key) + '\t' + str(count_both_dic[key]) + '\n'
			count_both_file.write(wboth)
	
		count_bc1_file = open(count_bc1_dir, 'w')
		for key1 in count_bc1_dic:
			wbc1 = str(key1) + '\t' + str(count_bc1_dic[key1]) + '\n'
			count_bc1_file.write(wbc1)
		
		count_bc2_file = open(count_bc2_dir, 'w')
		for key2 in count_bc2_dic:
			wbc2 = str(key2) + '\t' + str(count_bc2_dic[key2]) + '\n'
			count_bc2_file.write(wbc2)
		
		
		count_both_file.close()
		count_bc1_file.close()
		count_bc2_file.close()
			
		sorted_count_both_dic = sorted(count_both_dic.items(), key=lambda x: x[1], reverse=True)
		sorted_count_bc1_dic = sorted(count_bc1_dic.items(), key=lambda y: y[1], reverse=True)
		sorted_count_bc2_dic = sorted(count_bc2_dic.items(), key=lambda z: z[1], reverse=True)
	
	
		count_bc_file = open(count_bc_dir, 'w')
		both_count =  str.join('\t', [str(k)  for k  in sorted_count_both_dic])
		BOTH = 'Both_BC\t' + both_count + '\n'
		count_bc_file.write(BOTH)
	
		bc1_count =   str.join('\t', [str(k1)  for k1  in sorted_count_bc1_dic])
		bc1_c = 'BC1\t' + bc1_count + '\n'
		count_bc_file.write(bc1_c)
	
		bc2_count =   str.join('\t', [str(k2)  for k2  in sorted_count_bc2_dic])
		bc2_c = 'BC2\t' + bc2_count + '\n'
		count_bc_file.write(bc2_c)
	
		count_bc_file.close()
	
	
	
	
	if most_freq_BC == True:
	
	
		count_bc1_file = open(count_bc1_dir, 'r')
		count_bc2_file =  open(count_bc2_dir, 'r')
		
		previous_freq1 = 1
		previous_freq2 = 1
		h_bc1 = 'None'
		h_freq1 = 'NS'
		h_bc2 = 'None'
		h_freq2 = 'NS'
		
		count_line1 = 0
		for line1 in count_bc1_file:
			bc1, freq1 = line1.rstrip().split()
			count_line1 += 1
			if int(freq1) > (previous_freq1 + 1):
				h_bc1 = bc1
				h_freq1 = freq1
				previous_freq1 = int(freq1)
		if count_line1 ==1 and h_bc1 == 'None':
			h_bc1 = bc1
			h_freq1 = freq1
		
		count_line2 = 0
		for line2 in count_bc2_file:
			bc2, freq2 = line2.rstrip().split()
			count_line2
			if int(freq2) > (previous_freq2 + 1):
				h_bc2 = bc2
				h_freq2 = freq2
				previous_freq2 = int(freq2)
		if count_line2 ==1 and h_bc2 == 'None':
			h_bc2 = bc2
			h_freq2 = freq2
			
	
		hw = str.join('\t',(strain, h_bc2, h_bc1, h_freq2, h_freq1))
		h_freq_file.write(hw + '\n')
	
		count_bc1_file.close()
		count_bc2_file.close()
	
	
h_freq_file.close()
	
	
	
	
	
	
	
	
	
	
	
	
			