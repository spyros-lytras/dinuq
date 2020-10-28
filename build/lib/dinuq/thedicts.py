
####       DICTIONARIES        ####


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



############################



####       LISTS        ####    

#non-informative dinucleotide positions that will be excluded
noninfo = ['CpCpos1', 'CpApos1', 'GpCpos1', 'GpGpos1', 'GpUpos1', 'GpApos1', 'UpGpos1', 'UpApos1', 'ApCpos1', 'ApUpos1', 'ApApos1']    

stop = ["TAA", "TAG", "TGA"]
    
RNA_codons = ['UUU', 'UCU', 'UAU', 'UGU', 'UUC', 'UCC', 'UAC', 'UGC', 'UUA', 'UCA', 'UAA', 'UGA', 'UUG', 'UCG', 'UAG', 'UGG', 'CUU', 'CCU', 'CAU', 'CGU', 'CUC', 'CCC', 'CAC', 'CGC', 'CUA', 'CCA', 'CAA', 'CGA', 'CUG', 'CCG', 'CAG', 'CGG', 'AUU', 'ACU', 'AAU', 'AGU', 'AUC', 'ACC', 'AAC', 'AGC', 'AUA', 'ACA', 'AAA', 'AGA', 'AUG', 'ACG', 'AAG', 'AGG', 'GUU', 'GCU', 'GAU', 'GGU', 'GUC', 'GCC', 'GAC', 'GGC', 'GUA', 'GCA', 'GAA', 'GGA', 'GUG', 'GCG', 'GAG', 'GGG']

DNA_codons = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG', 'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']


############################