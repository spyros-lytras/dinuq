from Bio import SeqIO
from Bio.Seq import Seq
import random
import math

from dinuq.make_ntcont import eprime, ntcont
from dinuq.thedicts import syco, noninfo
from dinuq.checks import removeNs



###########################################################
###################    corrected SDU    ###################
###########################################################




#dinucl should be a list like: ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA']
#position should also be a list: ['pos1', 'pos2', 'bridge']
#boots can be any integer, recommend 100-1000
#custom_nt can be a dictionary of custom nucleotide compositions matching the sequence accessions (like output of ntcont())
    #{'acc1': {'G': fG, 'C': fC, 'A': fA, 'T': fT}, 'acc2': ...}

	
def SDUc(fasta_file, dinucl, position = ['pos1', 'pos2', 'bridge'], boots = 'none', custom_nt = 'none'):    


    ####    ALLOW CUSTOM NT COMPOSITION    ####
    
    #allow for custom nucleotide dictionary
    if custom_nt != 'none':
        
        ####    troubleshooting messages   ####
        for n in list(custom_nt):
            #check if all nts are there
            if set(custom_nt[n]) != {'A', 'G', 'T', 'C'}:
                print("\n\n\\(*'o'*)/\tOops! %s doesn't have the right nucleotides...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return
            #check if frequencies sum to 1
            if round(sum([custom_nt[n][d] for d in list(custom_nt[n])]), 10) != 1:
                print("\n\n\\(*'o'*)/\tOops! The %s nucleotide frequencies don't add up to 1...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return

        ntcontdic = custom_nt
        
    else:
        #calculate nucleotide composition dictionaries for all sequences
        ntcontdic = ntcont(fasta_file)

    ###########################################


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
        
        
        #remove Ns that might be in the sequence
        seq = removeNs(rec.id, seq)
        
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
        

        

        #the amino acid sequence list is independent of the dinucleotide CDS position
        aalist = []
                
        #for each amino acid in the sequence (range of the length of the sequence)
        for a in range(len(aa)):
        
            #append it to the aa list
            aalist.append(aa[a])       
        


        
        ####    BOOTS    ####


        #this will calculate the expected codon frequencies based on the single nucleotide composition
        #we will use this for weighing the synonymous codon sampling for the model sequence
            #position and dinucleotide options are given arbitrarily (null), since they are not used for returning the codon freqs
        codfreqs = eprime(ntcontdic[rec.id], 'null', 'null', extras = 'codonfreq')


        #list with actual sequence as the first item, followed by all the model sequences
        all_the_seqs = [seq]
        
        #if the boots option has been selected
        if boots != 'none': 
        
            #define how many random samples you take
            for i in range(boots):

                modelseq = "" 

                #for each amino acid in the translated sequence
                for a in aa:
                
                    #in order to account for single nucleotide composition in the sequence modelling
                    #we will use the expected codon frequencies as weights for 'random' sampling in a list
                    theweights = []            
                    for sc in syco[a]:
                        theweights.append(codfreqs[a][sc])
                    
                    #weighted 'random' selection of one possible synonymous codon for the amino acid
                    ran_cod = random.choices(syco[a], weights = theweights, k=1)
                                
                    #add the 'randomly' selected codon to the model sequence
                    modelseq = modelseq + ran_cod[0]
                    
                #add model sequence to the sequence list
                all_the_seqs.append(modelseq)
        

##########################


        #for each dinucleotide provided in the argument list
        for dinuc in dinucl:
        
            #for each position provided in the argument list
            for pos in position:
            
            
            
        ####        PARAMETERS        ####
        
                #define the dinucleotide without the linking p
                thedinucleotide = dinuc.replace('p', '')
                #don't forget to use DNA instead of RNA code!!
                thedinucleotide = thedinucleotide.replace('U', 'T')
                
                #calculate the corrected e values for this sequence composition,
                #dinucleotide and position using eprime
                thedictionary = eprime(ntcontdic[rec.id], dinuc, pos)
                
              

        ##########################


                        

    ####       ALL        ####
                      
                
                #bridge

                
                if pos == 'bridge':
                
                    #define the name of dinucleotide and position you are executing in this loop 
                    name = str(str(dinuc) + str(pos))
                    
                
                    #use this counter to separate first sequence in the list from the model ones
                    c = 0
                    #calculate wsdu for the observed sequence and all the model sequences
                    for sequ in all_the_seqs:
  
                    
                        #create lists for bridge dinucleotides
                        bdint = []                    

                        #define the last dinucleotide as the sequuence length -3
                        lastdi = len(sequ) - 3
                        
                        #for nucleotide starting at position 2 (3rd position of first codon) and ending at the last position with a step of 3 (one codon)
                        for d in range(2,lastdi, 3):
                        
                            #store dinucleotide
                            dint = str(sequ[d] + sequ[d+1])
                            #append it to the bridge dinucleotide list
                            bdint.append(dint)


                        #create dictionary that will get keys of amino acid pairs and values of lists of synonymous bridge dinucleotides for this pair
                        obsbridge = {}
                        
                        #for each item in the dn list
                        for i in range(len(bdint)):
                        
                            #create the amino acid pair
                            bridgeaa = str(aalist[i] + aalist[i+1])
                            
                            #if this amino acid pair has not been added to the list already, 
                            #update the dictionary with the pair as key and the relevant dn in the value list
                            if bridgeaa not in obsbridge:
                                    obsbridge.update({bridgeaa:[bdint[i]]})
                                    
                            #if it has already been added as a key, append the new dn to the list of that key
                            else:
                                    obsbridge[bridgeaa].append(bdint[i])
                       


                        #set all the variables to 0
                        wsdu = 0
                        k = 0

                        #for each amino acid pair key in the observed synonymous dn dictionary 
                        for dia in obsbridge:
                                
                            #if the pair is in the dictionary -> can use the dn of interest synonymously
                            if dia in thedictionary:

                                #store the synonymous dn list for that aa pair
                                dibridgelist = obsbridge[dia]
                                
                                #k is counting the total number of amino acid pairs of interest observed
                                k = k + len(dibridgelist)
                                
                                #count the number of dn of interest in the list 
                                num_obs = dibridgelist.count(thedinucleotide)
                         
                                #make the proportion of the dn of interest over all synonymous dns used
                                prop_obs = num_obs/len(dibridgelist)
                                
                                #store the expected proportion of the dn under expected synonymous usage for the aa pair
                                prop_exp = thedictionary[dia]
                                
                        ####    WEIGHTED    ####
                                
                                #store the wsdu for this aa pair as proportion observed over proportion expected dn timed by the number of occurences of this aa pair
                                wsdu_i = (prop_obs/prop_exp)*len(dibridgelist)
                                
                                #sum all wsdu values in the loop
                                wsdu = wsdu + wsdu_i
                        
                        #calculate the wsdu as the wsdu sum over the total cumulative number of aa pairs, instead of the wsdu values in the loop, k
                        # wsdu = wsdu/k
                        #takes care of division by 0 issue
                        if k != 0:
                            wsdu = wsdu/k
                        else:
                            wsdu = 'NA'

                        ########################
                        
                        
                        #if this is the first sequence in the seq list, add that outside the inner list (observed wsdu)
                        if c == 0:
                            if boots == 'none':
                                results.update({name:[wsdu]})
                            else:
                                results.update({name:[wsdu,[]]}) 
                                
                        #if this is not the first sequence, add it inside the inner list (modelled wsdu)
                        if c > 0:
                            results[name][1].append(wsdu)
                        
                        c = c + 1
                        
                        
                            
 


                   
                #pos1 or 2

                

                if (pos == 'pos1' or pos == 'pos2'):
                    
                    
                    #define pos_start variable from pos string
                    pos_start = int(pos.replace('pos', ''))-1
 
                    name = str(str(dinuc) + str(pos))
                    
                    if name not in noninfo:
                    
                        c = 0
                        
                        
                        for sequ in all_the_seqs:
                        
                            posdint = []
                            for d in range(pos_start,len(sequ), 3):
                                dint = str(sequ[d] + sequ[d+1])
                                posdint.append(dint)
                                       
                        
                            obspos = {}
                            for i in range(len(posdint)):
                                aa = aalist[i]
                                if aa not in obspos:
                                        obspos.update({aa:[posdint[i]]})
                                else:
                                        obspos[aa].append(posdint[i])
                                
                            
                            
                            wsdu = 0
                            k = 0
                            for aa in obspos:
                                if aa in thedictionary:
                                    di2list = obspos[aa]
                                    k = k + len(di2list)
                                    num_obs = di2list.count(thedinucleotide)
                                    prop_obs = num_obs/len(di2list)
                                    prop_exp = thedictionary[aa]
                                    wsdu_i = (prop_obs/prop_exp)*len(di2list)
                                    wsdu = wsdu + wsdu_i
                            
                            #takes care of division by 0 issue
                            if k != 0:
                                wsdu = wsdu/k
                            else:
                                wsdu = 'NA'
                            
                            if c == 0:
                                if boots == 'none':
                                    results.update({name:[wsdu]})
                                else:
                                    results.update({name:[wsdu,[]]})
                                
                            if c > 0:
                                results[name][1].append(wsdu)                            
                            c = c + 1
                                    
        all_results.update({rec.id:results})

        
    ##########################


    return all_results
    

##########################   









#dinucl should be a list like: ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA']
#position should also be a list: ['pos1', 'pos2', 'bridge']

	
def RSDUc(fasta_file, dinucl, position = ['pos1', 'pos2', 'bridge'], boots = 'none', custom_nt = 'none'):    


    ####    ALLOW CUSTOM NT COMPOSITION    ####
    
    #allow for custom nucleotide dictionary
    if custom_nt != 'none':
        
        ####    troubleshooting messages   ####
        for n in list(custom_nt):
            #check if all nts are there
            if set(custom_nt[n]) != {'A', 'G', 'T', 'C'}:
                print("\n\n\\(*'o'*)/\tOops! %s doesn't have the right nucleotides...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return
            #check if frequencies sum to 1
            if round(sum([custom_nt[n][d] for d in list(custom_nt[n])]), 10) != 1:
                print("\n\n\\(*'o'*)/\tOops! The %s nucleotide frequencies don't add up to 1...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return

        ntcontdic = custom_nt
        
    else:
        #calculate nucleotide composition dictionaries for all sequences
        ntcontdic = ntcont(fasta_file)

    ###########################################
    
    

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

        #remove Ns that might be in the sequence
        seq = removeNs(rec.id, seq)            
        
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
        

        

        #the amino acid sequence list is independent of the dinucleotide CDS position
        aalist = []
                
        #for each amino acid in the sequence (range of the length of the sequence)
        for a in range(len(aa)):
        
            #append it to the aa list
            aalist.append(aa[a])       
        


        
        ####    BOOTS    ####


        #this will calculate the expected codon frequencies based on the single nucleotide composition
        #we will use this for weighing the synonymous codon sampling for the model sequence
            #position and dinucleotide options are given arbitrarily (null), since they are not used for returning the codon freqs
        codfreqs = eprime(ntcontdic[rec.id], 'null', 'null', extras = 'codonfreq')


        #list with actual sequence as the first item, followed by all the model sequences
        all_the_seqs = [seq]
        
        #if the boots option has been selected
        if boots != 'none':
        
            #define how many random samples you take
            for i in range(boots):

                modelseq = "" 

                #for each amino acid in the translated sequence
                for a in aa:
                
                    #in order to account for single nucleotide composition in the sequence modelling
                    #we will use the expected codon frequencies as weights for 'random' sampling in a list
                    theweights = []            
                    for sc in syco[a]:
                        theweights.append(codfreqs[a][sc])
                    
                    #weighted 'random' selection of one possible synonymous codon for the amino acid
                    ran_cod = random.choices(syco[a], weights = theweights, k=1)
                                
                    #add the 'randomly' selected codon to the model sequence
                    modelseq = modelseq + ran_cod[0]
                    
                #add model sequence to the sequence list
                all_the_seqs.append(modelseq)
        

##########################


        #for each dinucleotide provided in the argument list
        for dinuc in dinucl:
        
            #for each position provided in the argument list
            for pos in position:
                
                
                
            
        ####        PARAMETERS        ####
        
                #define the dinucleotide without the linking p
                thedinucleotide = dinuc.replace('p', '')
                #don't forget to use DNA instead of RNA code!!
                thedinucleotide = thedinucleotide.replace('U', 'T')
                
                #calculate the corrected e values for this sequence composition,
                #dinucleotide and position using eprime
                thedictionary = eprime(ntcontdic[rec.id], dinuc, pos)
                
              

        ##########################


                        

    ####       ALL        ####


                #bridge

                
                if pos == 'bridge':
                
                    #define the name of dinucleotide and position you are executing in this loop 
                    name = str(str(dinuc) + str(pos))
                    
                    #making sure that at least one informative aa is in the sequence, so that the denominator is not 0
                    noz = 0
                    for theaa in list(thedictionary):
                        doubleaalist = [x+y for x,y in zip(aalist[0::2], aalist[1::2])]
                        if theaa in doubleaalist:
                            noz = noz + 1                    
                    
                    if noz >0:
                        
                        #use this counter to separate first sequence in the list from the model ones
                        c = 0
                        #calculate sdu for the observed sequence and all the model sequences
                        for sequ in all_the_seqs:
      
                        
                            #create lists for bridge dinucleotides
                            bdint = []                    

                            #define the last dinucleotide as the sequuence length -3
                            lastdi = len(sequ) - 3
                            
                            #for nucleotide starting at position 2 (3rd position of first codon) and ending at the last position with a step of 3 (one codon)
                            for d in range(2,lastdi, 3):
                            
                                #store dinucleotide
                                dint = str(sequ[d] + sequ[d+1])
                                #append it to the bridge dinucleotide list
                                bdint.append(dint)


                            #create dictionary that will get keys of amino acid pairs and values of lists of synonymous bridge dinucleotides for this pair
                            obsbridge = {}
                            
                            #for each item in the dn list
                            for i in range(len(bdint)):
                            
                                #create the amino acid pair
                                bridgeaa = str(aalist[i] + aalist[i+1])
                                
                                #if this amino acid pair has not been added to the list already, 
                                #update the dictionary with the pair as key and the relevant dn in the value list
                                if bridgeaa not in obsbridge:
                                        obsbridge.update({bridgeaa:[bdint[i]]})
                                        
                                #if it has already been added as a key, append the new dn to the list of that key
                                else:
                                        obsbridge[bridgeaa].append(bdint[i])
                           
         

                            #set all the variables to 0
                            rsdu = 0
                            sdu = 0
                            sdu_max = 0

                            #for each amino acid pair key in the observed synonymous dn dictionary 
                            for dia in obsbridge:
                                    
                                #if the pair is in the dictionary -> can use the dn of interest synonymously
                                if dia in thedictionary:
                                    
                                    #store the synonymous dn list for that aa pair
                                    dibridgelist = obsbridge[dia]
                                    
                                    #count the number of dn of interest in the list 
                                    num_obs = dibridgelist.count(thedinucleotide)
                             
                                    #make the proportion of the dn of interest over all synonymous dns used
                                    prop_obs = num_obs/len(dibridgelist)
                                    
                                    #store the expected proportion of the dn under expected synonymous usage for the aa pair
                                    prop_exp = thedictionary[dia]
                                    
                            ####    WEIGHTED    ####
                                    
                                    #store the sdu for this aa pair as proportion observed over proportion expected dn timed by the number of occurences of this aa pair
                                    sdu_i = (prop_obs/prop_exp)*len(dibridgelist)
                                    
                                    #define the rsdu denominator
                                    sdu_max_i = (1/prop_exp)*len(dibridgelist)
                                    
                                    #sum all sdu values in the loop
                                    sdu = sdu + sdu_i
                                    
                                    #sum all the denominators
                                    sdu_max = sdu_max + sdu_max_i

                            rsdu = sdu/sdu_max
                            

                            ########################
                            
                            
                            #if this is the first sequence in the seq list, add that outside the inner list (observed rsdu)
                            if c == 0:
                                if boots == 'none':
                                    results.update({name:[rsdu]})
                                else:
                                    results.update({name:[rsdu,[]]}) 
                            
                            #if this is not the first sequence, add it inside the inner list (modelled rsdu)
                            if c > 0:
                                results[name][1].append(rsdu)
                            
                            c = c + 1
                            
                    else:
                        if boots == 'none':
                            results.update({name:['NA']})
                        else:
                            results.update({name:['NA',['NA']]})                            
                            
                                
 


                   
                #pos1 or 2

                

                if (pos == 'pos1' or pos == 'pos2'):


                    #define pos_start variable from pos string
                    pos_start = int(pos.replace('pos', ''))-1
                                        
                    name = str(str(dinuc) + str(pos))
                    
                    if name not in noninfo:
                    
                        c = 0
                        
                        noz = 0
                        for thea in list(thedictionary):
                            if thea in aalist:
                                noz = noz + 1
                                
                        if noz > 0:
                        
                            for sequ in all_the_seqs:
                            
                                posdint = []
                                for d in range(pos_start,len(sequ), 3):
                                    dint = str(sequ[d] + sequ[d+1])
                                    posdint.append(dint)
                                           
                            
                                obspos = {}
                                for i in range(len(posdint)):
                                    aa = aalist[i]
                                    if aa not in obspos:
                                            obspos.update({aa:[posdint[i]]})
                                    else:
                                            obspos[aa].append(posdint[i])
                                    
                                
                                rsdu = 0
                                sdu = 0
                                sdu_max = 0
                                for aa in obspos:
                                    if aa in thedictionary:
                                        di2list = obspos[aa]
                                        num_obs = di2list.count(thedinucleotide)
                                        prop_obs = num_obs/len(di2list)
                                        prop_exp = thedictionary[aa]
                                        sdu_i = (prop_obs/prop_exp)*len(di2list)
                                        sdu_max_i = (1/prop_exp)*len(di2list)
                                        sdu = sdu + sdu_i
                                        sdu_max = sdu_max + sdu_max_i
                                rsdu = sdu/sdu_max
                                
                                if c == 0:
                                    if boots == 'none':
                                        results.update({name:[rsdu]})
                                    else:
                                        results.update({name:[rsdu,[]]}) 
                                if c > 0:
                                    results[name][1].append(rsdu)                            
                                c = c + 1
                        
                        else:
                            if boots == 'none':
                                results.update({name:['NA']})
                            else:
                                results.update({name:['NA',['NA']]})
                        
        all_results.update({rec.id:results})

    ##########################


    return all_results
    
   
##########################   








   
##########################   
   
   
def dict_to_tsv(dictionary, output_name, sep = '\t', error = 'none'):

    table_out = "acc" + sep

    first_accession = list(dictionary)[0]
    
    #with this if statement even if no error option is specified it will default to
    #95% CI if the dictionary includes uncertainty calculations
    if ( error == 'none' and len(dictionary[first_accession][list(dictionary[first_accession])[0]]) > 1 ):
        error = '95CI'
    
    #this takes care of the table's header
    for i in range(len(list(dictionary[first_accession]))):
    
        this_position = list(dictionary[first_accession])[i]
        
        if error == 'none':
        
            error_head = ''

        if error == 'stdev':
        
            error_head = str(this_position + '_lowSTDEV' + sep + this_position + '_highSTDEV' + sep)

        if error == '95CI':
        
            error_head = str(this_position + '_low95CI' + sep + this_position + '_high95CI' + sep)
        
        if error == 'extrema':
        
            error_head = str(this_position + '_errorMin' + sep + this_position + '_errorMax' + sep)
        
        add = str(this_position + sep + error_head)
        
        table_out = table_out + add
        
    table_out = str(table_out[:-1] + '\n')

    #for each accession key in the outer dictionary
    for i in range(len(list(dictionary))):
    
        acc = list(dictionary)[i]
        
        acc_dict = dictionary[acc]
        
        table_out = str(table_out + acc + sep)
        
        #for each position key in the inner dictionary
        for p in range(len(acc_dict)):
        
            pos_list = acc_dict[list(acc_dict)[p]]
            
            #the observed sequence's value 
            the_value = pos_list[0]
                       
            table_out = table_out + str(str(the_value) + sep)
            
            if error != 'none':
                
                #takes care of not available RSDU values (none of the informative values amino acids in the sequence)
                if the_value != 'NA':
                    
                    #list of bootstrap sequence calculation
                    model_values = pos_list[1]
                    
                    #sample size
                    n = len(model_values)
                    
                    #mean
                    m = (sum(model_values)/n)
     
                    #calculate standard deviation
                    sdnum = 0
                    for s in model_values:
                        thenum = (s-m)**2
                        sdnum = sdnum + thenum
                    stdev = math.sqrt(sdnum/(n-1))
                         
                    if error == 'stdev':

                        error_add = str(str(m-stdev) + sep + str(m+stdev) + sep)
                        
                        table_out = table_out + error_add
                    
                    if error == '95CI':
                        
                        error_add = str(str(m-1.96*stdev) + sep + str(m+1.96*stdev) + sep)
                        
                        table_out = table_out + error_add                
                    
                    if error == 'extrema':
                    
                        #maximum
                        maxi = max(model_values)
                        #minimum
                        mini = min(model_values)
                        
                        error_add = str(str(mini) + sep + str(maxi) + sep)
                        
                        table_out = table_out + error_add 
                
                else:
                
                    table_out = str(table_out + 'NA' + sep + 'NA' + sep)
                
        table_out = str(table_out[:-1] + '\n')
        
    table_out = table_out[:-1]

    out = open(output_name, "w+")
    out.write(table_out)
    out.close()    
    
    return print(table_out)


