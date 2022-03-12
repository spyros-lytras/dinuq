from Bio import SeqIO
from Bio.Seq import Seq


def cod_check(fasta_file, errors = 'critical'):
    
	#parse the multiple sequence fasta file
    records = SeqIO.parse(fasta_file, "fasta")	

    problems = []


    #for each sequence
    for rec in records:

        #save sequence as string 
        seq = str(rec.seq)
                
        #ungap by simply removing dashes and tilts from the string
        seq = seq.replace('-', '')
        seq = seq.replace('~', '')
        
        #make everything upper case from the start
        seq = seq.upper()
        
        
        #make sure it's a coding sequence
        if len(seq)%3 != 0:
            print("\n\\(*'o'*)/\tOops! Sequence %s has length not a multiple of 3...\n"%str(rec.id))
            
            problems.append(rec.id)
            
        else:
                
            #save amino acid sequence as string
            aa = str(Seq(seq).translate())
            
            #remove stop codons at the end
            if aa[-1] == '*':
                aa = aa[:-1]
                seq = seq[:-3]
            
            #check for internal stop codons
            if '*' in aa:
                print("\n\\(*'o'*)/\tOops! Sequence %s has internal stop codons...\n"%str(rec.id))
                                
                problems.append(rec.id)
                
            
        #check for ambiguous bases
            #non-critical
        if (set(seq) != {'C', 'T', 'A', 'G'}) and (errors == 'all'):
            print("\n\\(*'o'*)/\tOops! Sequence %s has ambiguous nucleotides...\n"%str(rec.id))
            
            problems.append(rec.id)                
    
    return problems



#this is an internal function removing codons with ambiguous bases from a sequence
    #!! N removal could skew the bridge dinucleotide estimate !!
        #comes with a warning message
def removeNs(theid, seq):
    
    noNseq = ''

    if set(seq) != {'C', 'T', 'A', 'G'}:
        
        #warning
        print("\n\\(*'o'*)/\tSequence %s has ambiguous nucleotides!\n\tAmbiguous codons will be removed, but that might affect bridge dinucleotide calculations."%str(theid))
                
        for c in range(0,len(seq), 3):
            
            cod = seq[c:c+3]
            if set(cod) <= {'C', 'T', 'A', 'G'}:
                noNseq = noNseq + cod
                
    else:
        
        noNseq = seq
        
    return noNseq
                
                
    