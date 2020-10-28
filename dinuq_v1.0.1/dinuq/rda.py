from Bio import SeqIO
from Bio.Seq import Seq


###########################################################
###################         RDA         ###################
###########################################################





#non-informative dinucleotide positions that will be excluded
noninfo = ['CpCpos1', 'CpApos1', 'GpCpos1', 'GpGpos1', 'GpUpos1', 'GpApos1', 'UpGpos1', 'UpApos1', 'ApCpos1', 'ApUpos1', 'ApApos1']    


#dinucl should be a list like: ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA']
#position should also be a list: ['pos1', 'pos2', 'bridge', 'all']
	
def RDA(fasta_file, dinucl, position = ['all']):    



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
        

        #the amino acid sequence list is independent of the dinucleotide CDS position
        aalist = []
                
        #for each amino acid in the sequence (range of the length of the sequence)
        for a in range(len(aa)):
        
            #append it to the aa list
            aalist.append(aa[a])       
        
        
        
        ##########################


                

        ####   CALCULATIONS    ####  


        #calculate one RDA value for the entire sequence
        
        if 'all' in position:
        
           
            #for each dinucleotide provided in the argument list
            for dinuc in dinucl:
            
           
                
        ####        PARAMETERS        ####
        
        #define the parameters mapping to each dinucleotide and position
        
                if dinuc == 'CpC':
                    thedinucleotide = 'CC'

                if dinuc == 'CpG':
                    thedinucleotide = 'CG'

                if dinuc == 'CpU':
                    thedinucleotide = 'CT'

                if dinuc == 'CpA':
                    thedinucleotide = 'CA'

                if dinuc == 'GpC':
                    thedinucleotide = 'GC'

                if dinuc == 'GpG':
                    thedinucleotide = 'GG'

                if dinuc == 'GpU':
                    thedinucleotide = 'GT'

                if dinuc == 'GpA':
                    thedinucleotide = 'GA'

                if dinuc == 'UpC':
                    thedinucleotide = 'TC'

                if dinuc == 'UpG':
                    thedinucleotide = 'TG'

                if dinuc == 'UpU':
                    thedinucleotide = 'TT'

                if dinuc == 'UpA':
                    thedinucleotide = 'TA'

                if dinuc == 'ApC':
                    thedinucleotide = 'AC'

                if dinuc == 'ApG':
                    thedinucleotide = 'AG'

                if dinuc == 'ApU':
                    thedinucleotide = 'AT'

                if dinuc == 'ApA':
                    thedinucleotide = 'AA'


            ##########################                    
                                             

                #create lists for bridge dinucleotides
                dint = []                    

                #define the last dinucleotide as the sequence length
                lastdi = len(seq) - 1
                
                #for nucleotide starting at position 2 (3rd position of first codon) and ending at the last position with a step of 3 (one codon)
                for d in range(0,lastdi):
                
                    #store dinucleotide
                    di = str(seq[d] + seq[d+1])
                    #append it to the bridge dinucleotide list
                    dint.append(di)


                thedinucleotideone = thedinucleotide[0]
                thedinucleotidetwo = thedinucleotide[1]

                #calculate frequencies
                freq_one = seq.count(thedinucleotideone)/len(seq)
                freq_two = seq.count(thedinucleotidetwo)/len(seq)
                freq_all = dint.count(thedinucleotide)/len(dint)
                
                #calculate rda
                rda = freq_all/(freq_one*freq_two)

                
                name = str(dinuc)
                
                results.update({name:rda})



            all_results.update({rec.id:results})
                
                
                

        ####       POSITIONS        ####
                
        
        #calculate a separate RDA value for each frame position
        else:    
           

            #for each dinucleotide provided in the argument list
            for dinuc in dinucl:
            
           
                
        ####        PARAMETERS        ####
        
        #define the parameters mapping to each dinucleotide and position
        
                if dinuc == 'CpC':
                    thedinucleotide = 'CC'

                if dinuc == 'CpG':
                    thedinucleotide = 'CG'

                if dinuc == 'CpU':
                    thedinucleotide = 'CT'

                if dinuc == 'CpA':
                    thedinucleotide = 'CA'

                if dinuc == 'GpC':
                    thedinucleotide = 'GC'

                if dinuc == 'GpG':
                    thedinucleotide = 'GG'

                if dinuc == 'GpU':
                    thedinucleotide = 'GT'

                if dinuc == 'GpA':
                    thedinucleotide = 'GA'

                if dinuc == 'UpC':
                    thedinucleotide = 'TC'

                if dinuc == 'UpG':
                    thedinucleotide = 'TG'

                if dinuc == 'UpU':
                    thedinucleotide = 'TT'

                if dinuc == 'UpA':
                    thedinucleotide = 'TA'

                if dinuc == 'ApC':
                    thedinucleotide = 'AC'

                if dinuc == 'ApG':
                    thedinucleotide = 'AG'

                if dinuc == 'ApU':
                    thedinucleotide = 'AT'

                if dinuc == 'ApA':
                    thedinucleotide = 'AA'


            ##########################              

                  
                #for each position provided in the argument list
                for pos in position:
            
                #bridge

                    
                    if pos == 'bridge':
      
                        
                        #create lists for bridge dinucleotides
                        bdint = []                    

                        #define the last dinucleotide as the sequence length -3
                        lastdi = len(seq) - 3
                        
                        #for nucleotide starting at position 2 (3rd position of first codon) and ending at the last position with a step of 3 (one codon)
                        for d in range(2,lastdi, 3):
                        
                            #store dinucleotide
                            dint = str(seq[d] + seq[d+1])
                            #append it to the bridge dinucleotide list
                            bdint.append(dint)

                        
                        dbstr = ""
                        for d in bdint:
                            dbstr = dbstr + d   

                        thedinucleotideone = thedinucleotide[0]
                        thedinucleotidetwo = thedinucleotide[1]

                        freq_one = dbstr.count(thedinucleotideone)/len(dbstr)
                        freq_two = dbstr.count(thedinucleotidetwo)/len(dbstr)
                        freq_all = bdint.count(thedinucleotide)/len(bdint)
                        rda = freq_all/(freq_one*freq_two)

                        
                        name = str(str(dinuc) + str(pos))
                        
                        results.update({name:rda})
     


                       
                    #pos1 or 2

                    

                    if (pos == 'pos1' or pos == 'pos2'):
                        
                        name = str(str(dinuc) + str(pos))
                        
                        if name not in noninfo:
                        
                            if pos == 'pos1':
                                pos_start = 0
                            if pos == 'pos2':
                                pos_start = 1
                                
                            posdint = []
                            for d in range(pos_start,len(seq), 3):
                                dint = str(seq[d] + seq[d+1])
                                posdint.append(dint)
                                       
                            
                            dbstr = ""
                            for d in posdint:
                                dbstr = dbstr + d   

                            thedinucleotideone = thedinucleotide[0]
                            thedinucleotidetwo = thedinucleotide[1]

                            freq_one = dbstr.count(thedinucleotideone)/len(dbstr)
                            freq_two = dbstr.count(thedinucleotidetwo)/len(dbstr)
                            freq_all = posdint.count(thedinucleotide)/len(posdint)
                            rda = freq_all/(freq_one*freq_two)
                                                  
                            
                            results.update({name:rda})

                     
            
        all_results.update({rec.id:results})

    ##########################


    return all_results
    











####    TABLE   ####
   
   
   
def RDA_to_tsv(rda_dic, output_name):

    table_out = "acc\t"

    for i in range(len(list(rda_dic[list(rda_dic)[0]]))):
        add = str(list(rda_dic[list(rda_dic)[0]])[i] + '\t')
        table_out = table_out + add
        
    table_out = str(table_out + '\n')

    for i in range(len(list(rda_dic))):
        acc = str(list(rda_dic)[i] + '\t')
        dict = rda_dic[list(rda_dic)[i]]
        table_out = table_out + acc
        for r in range(len(dict)):
            this_number = str(str(dict[list(dict)[r]]) + '\t')
            table_out = table_out + this_number
        table_out = str(table_out + '\n')
        
    table_out = table_out[:-1]
    
    out = open(output_name, "w+")
    out.write(table_out)
    out.close()    
    
    return print(table_out) 
    
    
##########################    
