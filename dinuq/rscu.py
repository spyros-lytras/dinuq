from Bio import SeqIO
from Bio.Seq import Seq

from dinuq.thedicts import syco, cod_aa_dic, stop, DNA_codons, RNA_codons

###########################################################
###################         RSCU         ##################
###########################################################


	
def RSCU(fasta_file):    



    #define all results dictionary    
    all_results = {}



    #parse the multiple sequence fasta file
    records = SeqIO.parse(fasta_file, "fasta")	


    ####       FASTA        ####


    #for each sequence
    for rec in records:
    
        results = {}
        
        

        #save sequence as string 
        seq = str(rec.seq)
                
        #ungap by simply removing dashes and tilts from the string
        seq = seq.replace('-', '')
        seq = seq.replace('~', '')
        
        #make everything upper case from the start
        seq = seq.upper()
        
        
        #make sure it's a coding sequence
        if len(seq)%3 != 0:
            print("\n\n\\(*'o'*)/\tOops! Sequence %s has length not a multiple of 3...\n\n"%str(rec.id))
            return
        
        
        #save amino acid sequence as string
        aa = str(Seq(seq).translate())
        
        #remove stop codons at the end
        if aa[-1] == '*':
            aa = aa[:-1]
            seq = seq[:-3]
        
        #check for internal stop codons
        if '*' in aa:
            print("\n\n\\(*'o'*)/\tOops! Sequence %s has internal stop codons...\n\n"%str(rec.id))
            return



        
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
            if c not in stop:
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


def RSCU_to_tsv(RSCU_dic, output_name, sep = '\t'):

    table_out = ""

    for acc in list(RSCU_dic):
    
        #add the accession
        table_out = str(table_out + str(acc) + '\n')
        
        #need to add a return every 4 AAs, so add a counter
        rc = 1
        
        for i in range(len(DNA_codons)):
        
            codon = DNA_codons[i]
            
            if codon not in stop:
            
                value = RSCU_dic[acc][codon]
                
                if rc%4 != 0:
                    
                    table_out = str(table_out + cod_aa_dic[codon] + sep + RNA_codons[i] + sep + str(value) + sep)
                    
                else:
                
                    table_out = str(table_out + cod_aa_dic[codon] + sep + RNA_codons[i] + sep + str(value) + '\n')
                    
            else: 
                
                value = 'NA'
                
                if rc%4 != 0:
                    
                    table_out = str(table_out + cod_aa_dic[codon] + sep + RNA_codons[i] + sep + str(value) + sep)
                    
                else:
                
                    table_out = str(table_out + cod_aa_dic[codon] + sep + RNA_codons[i] + sep + str(value) + '\n')

            rc = rc + 1
            
       
        table_out = str(table_out + '\n\n')
        
        
    table_out = table_out[:-2]
    
    out = open(output_name, "w+")
    out.write(table_out)
    out.close()
    
    return print(table_out)
    
    


##########################   

    
    
