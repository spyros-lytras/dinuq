from Bio import SeqIO

from dinuq.thedicts import syco, stop




def ntcont(fasta):
    
    """parses sequences and makes a dictionary with
    single nt compositions for each seq in a fasta file"""

    ntcontdic = {}
            
    records = SeqIO.parse(fasta, "fasta")	


    ####       FASTA        ####


    #for each sequence
    for rec in records:

        # #ungap in case of an alignment
        # recungap = rec.seq.ungap("-")
        
        #save sequence as string 
        seq = str(rec.seq)
                
        #ungap by simply removing dashes and tilts from the string
        seq = seq.replace('-', '')
        seq = seq.replace('~', '')
        
        #make everything upper case from the start
        seq = seq.upper()
        
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

        ntcontdic.update({rec.id:{'G':Gcon, 'C':Ccon, 'A':Acon, 'T':Tcon}})
        
    return ntcontdic
    
    
    
    





def eprime(composition, dinucl, position, extras = 'eprime'):
    
    """
    this function calculates 3 different sets of values,
    all based on the nucleotide composition provided:
        
        i)stop codon correction factor
        
        ii)dict of corrected expected synonymous codon frequencies

        iii)dict of corrected expected dinucleotide frequencies for a
            specified dinucleotide and codon position
            
    this function is not standalone, but to be used by sduc and rsduc
    or used for multiple sequences by eprimeall
    """
    
    
    ntdic = composition
    
    
    ####    correction factor   #####
    #predefine the correction factor based on the stop codon frequencies
    
    #this will store the sum of all stop codon exp freqs
    fstop = 0
    
    for c in stop:
        
        fsa = ntdic[c[0]]

        fsb = ntdic[c[1]]

        fsc = ntdic[c[2]]

        fsco = fsa*fsb*fsc
        
        fstop = fstop + fsco
        
    #eq for calculating correction factor from stop codon frequencies (see OneNote)
    corr = 1/(1-fstop)
    
    
    ####    EXTRAS   ####
    #only return the correction factor
    if extras == 'correction':        
        return corr
    #####################
    
    
    #this will be the expected codon frequency dictionary
    sycofreq = {}
    
    #for each amino acid
    for a in syco:
    
        #define the new dict
        afreqs = {}
        
        #for each synonymous codon of the aa
        for c in syco[a]:
        
            #pos1 nucleotide freq
            fa = ntdic[c[0]]

            #pos2 nucleotide freq
            fb = ntdic[c[1]]
            
            #pos3 nucleotide freq            
            fc = ntdic[c[2]]
            
            #multiply to get exp freq of codon
            fco = fa*fb*fc
            
            #correct with the precalculated correction factor
            fcoc = fco*corr
            
            #update the dict with codon:corrected exp frequency
            afreqs.update({c:fcoc})
            
        #update the final dict with amino acid : {synonymous codon : corr exp freq, ...}
        sycofreq.update({a:afreqs})



    ####    EXTRAS   ####
    #only return the corrected expected synonymous codon frequencies
    if extras == 'codonfreq':
        return sycofreq
    #####################
    
    
                
    #Now! off to calculating the expected dinucleotide frequencies!!
    
    
    
    thedinucleotide = dinucl.replace('p', '')
    thedinucleotide = thedinucleotide.replace('U', 'T')

    pos_definitions = {'pos1':[0,2], 'pos2':[1,3], 'bridge':[]}

    pos = pos_definitions[position]


    
    #instead of precalculating everything define a single expect dict
    expect = {}
    
    #position definition step
    if position in ['pos1', 'pos2']:

        #for each amino acid 
        for a in list(sycofreq):
            
            #store the exp syn cod freq dict for it
            theaa = sycofreq[a]
            
            #make a list of the synonymous codons
            thesyncods = list(theaa)
            
            #list comprehension woohoo! 
            #stores all the synonymous codons for this aa that have the given dinucleotide in the given positions
            checker = [x for x in thesyncods if x[pos[0]:pos[1]] == thedinucleotide]
            
            #if the checker isn't empty, i.e. if the aa has syncods with the given parameters, 
            if len(checker) != 0: 
                
                #calculates the sum of frequencies of all synonymous codons   
                freqall = sum([theaa[c] for c in thesyncods])
                #calculates the sum of frequencies of the codons in checker
                freqdin = sum([theaa[c] for c in checker])
                #takes the ratio of the two, which equals to the expected synonymous dinucleotide frequency for this aa
                finfreq = freqdin/freqall
                
                #checks that aas are not added more than once
                if a not in list(expect):
                    #updates the expect dict with amino acid : expected synonymous dinucleotide frequency for this aa
                    expect.update({a:finfreq})
                
               
    
    #bridge
    if position == 'bridge':
        
        #a is first n is second
        #for each amino acid 
        for aa in list(sycofreq):
            
            #store the exp syn cod freq dict for it
            theaa = sycofreq[aa]
                
            #make a list of the synonymous codons
            thesyncodsa = list(theaa)
                
            
            #for the second amino acid
            for an in list(sycofreq):
            
                #store the exp syn cod freq dict for it
                thean = sycofreq[an]
                
                #make a list of the synonymous codons
                thesyncodsn = list(thean)
                
                #store the dipeptide
                dipep = aa + an            
                
                
                #this will check if the dinucleotide of interest can be given by the amino acid pair
                checker = [str(b[2] + m[0]) for b in thesyncodsa for m in thesyncodsn if str(b[2] + m[0]) == thedinucleotide]
                
                
                
                
                #if the checker isn't empty, i.e. if the aa has syncods with the given parameters, 
                if len(checker) != 0: 
                    
                    checkera = [x for x in thesyncodsa if x[2] == thedinucleotide[0]]
                    
                    checkern = [x for x in thesyncodsn if x[0] == thedinucleotide[1]]
                    
                    

                    #a
                    #calculates the sum of frequencies of all synonymous codons   
                    freqalla = sum([theaa[c] for c in thesyncodsa])
                    #calculates the sum of frequencies of the codons in checkera
                    freqdina = sum([theaa[c] for c in checkera])
                    #takes the ratio of the two, which equals to the expected synonymous dinucleotide frequency for this aa
                    finfreqa = freqdina/freqalla
                    
                    
                    #n
                    #calculates the sum of frequencies of all synonymous codons   
                    freqalln = sum([thean[c] for c in thesyncodsn])
                    #calculates the sum of frequencies of the codons in checkern
                    freqdinn = sum([thean[c] for c in checkern])
                    #takes the ratio of the two, which equals to the expected synonymous dinucleotide frequency for this an
                    finfreqn = freqdinn/freqalln
                    
                    
                    #the final bridge frequency will be the product of the frequency of synonymous codons for the first amino acid
                    #that have the first nt of the dinucleotide in their third codon position and the frequency of syn cods for the
                    #second aa that have the second nt of the dint in their first codon position
                    finfreq = finfreqa * finfreqn
           
                    
                    
                    #checks that aa pair are not added more than once
                    if dipep not in list(expect):
                        #updates the expect dict with amino acid : expected synonymous dinucleotide frequency for this aa pair
                        expect.update({dipep:finfreq})
    
    ####    EXTRAS    ####
    #this should be the default
    #returns the eprime values to be used for calculating sduc
    if extras == 'eprime':
        return expect
    ######################
    
    










def eprimeall(records, dinucl, position, extras = 'eprime'):
    
    """
    this function is a standalone version of eprime.
    It calculates the values below (chosen in extras option):
        
        i)stop codon correction factor
        
        ii)dict of corrected expected synonymous codon frequencies

        iii)dict of corrected expected dinucleotide frequencies for a
            specified dinucleotide and codon position
            
    given either:
        
        i) a fasta file with any number of sequences (can be aligned)
        
        ii) or a dictionary of custom nucleotide  composition dictionaries
    """
    
    
    #allow for custom nucleotide dictionary
    if type(records) is dict:
        
        
        ####    troubleshooting messages   ####
        for n in list(records):
            #check if all nts are there
            if set(records[n]) != {'A', 'G', 'T', 'C'}:
                print("\\(*'o'*)/\tOops! %s doesn't have the right nucleotides...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return
            #check if frequencies sum to 1
            if round(sum([records[n][d] for d in list(records[n])]), 10) != 1:
                print("\\(*'o'*)/\tOops! The %s nucleotide frequencies don't add up to 1...\n\nThe custom composition dictionary should be given as follows:\n\n\t{name:{'G':fG, 'C':fC, 'A':fA, 'T':fT},...}\n\n"%str(n))
                return        
        
        
        ntcontdic = records
        
    else:
        ntcontdic = ntcont(records)
    
    alldic = {}
    
    
    for rec in list(ntcontdic):
        
        comp = ntcontdic[rec]
                
        alldic.update({rec:eprime(comp, dinucl, position, extras)})
        
    return alldic

    