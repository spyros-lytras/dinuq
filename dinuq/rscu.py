from Bio import SeqIO
from Bio.Seq import Seq

###########################################################
###################         RSCU         ##################
###########################################################



###       DICTIONARIES        ####

syco = { 
    "C": ["TGT", "TGC"], 
    "D": ["GAT", "GAC"], 
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"], 
    "Q": ["CAA", "CAG"], 
    "M": ["ATG"], 
    "N": ["AAC", "AAT"], 
    "P": ["CCT", "CCG", "CCA", "CCC"], 
    "K": ["AAG", "AAA"], 
    "T": ["ACC", "ACA", "ACG", "ACT"], 
    "F": ["TTT", "TTC"], 
    "A": ["GCA", "GCC", "GCG", "GCT"], 
    "G": ["GGT", "GGG", "GGA", "GGC"], 
    "I": ["ATC", "ATA", "ATT"], 
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"], 
    "H": ["CAT", "CAC"], 
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"], 
    "W": ["TGG"], 
    "V": ["GTA", "GTC", "GTG", "GTT"], 
    "E": ["GAG", "GAA"], 
    "Y": ["TAT", "TAC"]} 


cod_aa_dic = {'TAA':'STOP', 'TAG': 'STOP', 'TGA': 'STOP', 'TGT': 'C', 'TGC': 'C', 'GAT': 'D', 'GAC': 'D', 'TCT': 'S', 'TCG': 'S', 'TCA': 'S', 'TCC': 'S', 'AGC': 'S', 'AGT': 'S', 'CAA': 'Q', 'CAG': 'Q', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N', 'CCT': 'P', 'CCG': 'P', 'CCA': 'P', 'CCC': 'P', 'AAG': 'K', 'AAA': 'K', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACT': 'T', 'TTT': 'F', 'TTC': 'F', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGT': 'G', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'ATC': 'I', 'ATA': 'I', 'ATT': 'I', 'TTA': 'L', 'TTG': 'L', 'CTC': 'L', 'CTT': 'L', 'CTG': 'L', 'CTA': 'L', 'CAT': 'H', 'CAC': 'H', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGG': 'R', 'AGA': 'R', 'TGG': 'W', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GAG': 'E', 'GAA': 'E', 'TAT': 'Y', 'TAC': 'Y'}

    
##################################
    


####       LISTS        ####    

    
stop_codons = ['TAA', 'TAG', 'TGA']

RNA_codons = ['UUU', 'UCU', 'UAU', 'UGU', 'UUC', 'UCC', 'UAC', 'UGC', 'UUA', 'UCA', 'UAA', 'UGA', 'UUG', 'UCG', 'UAG', 'UGG', 'CUU', 'CCU', 'CAU', 'CGU', 'CUC', 'CCC', 'CAC', 'CGC', 'CUA', 'CCA', 'CAA', 'CGA', 'CUG', 'CCG', 'CAG', 'CGG', 'AUU', 'ACU', 'AAU', 'AGU', 'AUC', 'ACC', 'AAC', 'AGC', 'AUA', 'ACA', 'AAA', 'AGA', 'AUG', 'ACG', 'AAG', 'AGG', 'GUU', 'GCU', 'GAU', 'GGU', 'GUC', 'GCC', 'GAC', 'GGC', 'GUA', 'GCA', 'GAA', 'GGA', 'GUG', 'GCG', 'GAG', 'GGG']

DNA_codons = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG', 'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']


############################



	
def RSCU(fasta_file):    



    #define all results dictionary    
    all_results = {}



    #parse the multiple sequence fasta file
    records = SeqIO.parse(fasta_file, "fasta")	


    ####       FASTA        ####


    #for each sequence
    for rec in records:
    
        results = {}


        #ungap in case of an alignment
        recungap = rec.seq.ungap("-")
        
        #save sequence as string 
        seq = str(recungap)
        recungap = rec.seq.ungap("~")
        
        #save sequence as string 
        seq = str(recungap)
        
        #make all uppercase
        seq = seq.upper()
        
        #make sure it's a coding sequence
        if len(seq)%3 != 0:
            print(str('\n\nSequence ' + rec.id + ' has length not a multiple of 3...\n\n'))
        
        
        #save amino acid sequence as string
        aa = str(recungap.translate())
        
        #remove stop codons at the end
        if aa[-1] == '*':
            aa = aa[:-1]
            seq = seq[:-3]
        
        #check for internal stop codons
        if '*' in aa:
            print(str('\n\nSequence ' + rec.id + ' has internal stop codons...\n\n'))
        
        #make an ordered list of all codons in the sequence
        cod = []
        for c in range(0,len(seq),3):
            codon = str(seq[c]+seq[c+1]+seq[c+2])
            cod.append(codon)        
        
        
        #make a dictionary with aa as keys and all observed synonymous codons as values
        cu_dic = {}

        for i in range(len(aa)):
            if aa[i] not in list(cu_dic):
                    cu_dic.update({aa[i]:[cod[i]]})
            else:
                    cu_dic[aa[i]].append(cod[i])
                    

        #calculate the rscu for each codon in the sequence        
        for c in list(cod_aa_dic):
            a = cod_aa_dic[c]
            if c not in stop_codons:
                if a in list(cu_dic):
                    x = cu_dic[a].count(c)
                    n = len(syco[a]) 
                    s_x = len(cu_dic[a])
                    the_rscu = (x/((1/n)*s_x))
                    results.update({c:the_rscu})
                else:
                    results.update({c:'NA'})
            
                       
    
        all_results.update({rec.id:results})

    ##########################


    return all_results
    
   
##########################   





####    TABLE   ####


def RSCU_to_tsv(RSCU_dic, output_name):

    table_out = ""

    for acc in list(RSCU_dic):
    
        #add the accession
        table_out = str(table_out + str(acc) + '\n')
        
        #need to add a return every 4 AAs, so add a counter
        rc = 1
        
        for i in range(len(DNA_codons)):
        
            codon = DNA_codons[i]
            
            if codon not in stop_codons:
            
                value = RSCU_dic[acc][codon]
                
                if rc%4 != 0:
                    
                    table_out = str(table_out + cod_aa_dic[codon] + '\t' + RNA_codons[i] + '\t' + str(value) + '\t')
                    
                else:
                
                    table_out = str(table_out + cod_aa_dic[codon] + '\t' + RNA_codons[i] + '\t' + str(value) + '\n')
                    
            else: 
                
                value = 'NA'
                
                if rc%4 != 0:
                    
                    table_out = str(table_out + cod_aa_dic[codon] + '\t' + RNA_codons[i] + '\t' + str(value) + '\t')
                    
                else:
                
                    table_out = str(table_out + cod_aa_dic[codon] + '\t' + RNA_codons[i] + '\t' + str(value) + '\n')

            rc = rc + 1
            
       
        table_out = str(table_out + '\n\n')
        
        
    table_out = table_out[:-2]
    
    out = open(output_name, "w+")
    out.write(table_out)
    out.close()
    
    return print(table_out)
    
    


##########################   

    
    
