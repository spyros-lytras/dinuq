from Bio import SeqIO

fasta_file = 'nCoVbat.pan6.fas'

ntcont = {}
		
records = SeqIO.parse(fasta_file, "fasta")	


####       FASTA        ####


#for each sequence
for rec in records:

	#ungap in case of an alignment
	recungap = rec.seq.ungap("-")
	
	#save sequence as string 
	seq = str(recungap)

	Gs = seq.count('G')
	Cs = seq.count('C')
	As = seq.count('A')
	Ts = seq.count('T')
	GCs = Gs + Cs
	GCcon = GCs/len(seq)
	Gcon = Gs/len(seq)
	Ccon = Cs/len(seq)
	Acon = As/len(seq)
	Tcon = Ts/len(seq)

	ntcont.update({rec.id:{'G':Gcon, 'C':Ccon, 'A':Acon, 'T':Tcon}})
	
print(ntcont)
