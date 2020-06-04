from Bio import SeqIO
from Bio.Seq import Seq
import random
import math

#total number of nucleotides used in the estimator
ntsamples = 19200000



syco_orig = { 
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



#### GC content ####

full_ntcont = {'RpShaanxi2011|Bat-R_pusillus|Shaanxi|JX993987|2011-09': {'G': 0.21082621082621084, 'C': 0.20509428842762176, 'A': 0.2843915343915344, 'T': 0.299687966354633}, 'HuB2013|Bat-Bat-R_sinicus|Hubei|KJ473814|2013-04': {'G': 0.21019623710297392, 'C': 0.2041270483512037, 'A': 0.28447636388158337, 'T': 0.301200350664239}, '279_2005|Bat-R_macrotis|Hubei|DQ648857|2004-11': {'G': 0.21014760767963417, 'C': 0.2010019837934165, 'A': 0.2851618977169564, 'T': 0.30368851080999293}, 'Rm1|Bat-R_macrotis|Hubei|DQ412043|2004-11': {'G': 0.20982217889676963, 'C': 0.20131769135096977, 'A': 0.28569027530337154, 'T': 0.30316985444888905}, 'JL2012|Bat-R_ferrumequinum|Jilin|KJ473811|2012-10': {'G': 0.20938802217859973, 'C': 0.20060612322209595, 'A': 0.284258015635224, 'T': 0.3057478389640803}, 'JTMC15|Bat-R_ferrumequinum|Jilin|KU182964|2013-10': {'G': 0.20896352699836584, 'C': 0.20089704808594971, 'A': 0.2831959945759883, 'T': 0.3069434303396961}, 'HeB2013|Bat-R_ferrumequinum|Hebei|KJ473812|2013-04': {'G': 0.2085724960092382, 'C': 0.20127025099344495, 'A': 0.284991339197772, 'T': 0.3051659137995449}, 'SX2013|Bat-R_ferrumequinum|Shanxi|KJ473813|2013-11': {'G': 0.20851294932283357, 'C': 0.20121516581242999, 'A': 0.28491904551780317, 'T': 0.30535283934693325}, 'Jiyuan-84|Bat-R_ferrumequinum|Henan-Jiyuan|KY770860|2012': {'G': 0.20845839937944757, 'C': 0.2020167953863276, 'A': 0.28461097433476107, 'T': 0.3049138308994638}, 'Rf1|Bat-R_ferrumequinum|Hubei-Yichang|DQ412042|2004-11': {'G': 0.20892658790265575, 'C': 0.20205998182368978, 'A': 0.2845602342724427, 'T': 0.30445319600121173}, 'GX2013|Bat-R_sinicus|Guangxi|KJ473815|2012-11': {'G': 0.2086005281026028, 'C': 0.2000274339014437, 'A': 0.28606700730427626, 'T': 0.30530503069167725}, 'Rp3|Bat-R_pearsoni|Guangxi-Nanning|DQ071615|2004-12': {'G': 0.2085687382297552, 'C': 0.19995964487489912, 'A': 0.2874293785310734, 'T': 0.30404223836427224}, 'Rf4092|Bat-R_ferrumequinum|Yunnan-Kunming|KY417145|2012-09-18': {'G': 0.21013126893301917, 'C': 0.20265903736115787, 'A': 0.2838774823291821, 'T': 0.30333221137664085}, 'Rs4231|Bat-R_sinicus|Yunnan-Kunming|KY417146|2013-04-17': {'G': 0.2086495198442012, 'C': 0.20085957961184608, 'A': 0.284668591766839, 'T': 0.3058223087771137}, 'WIV16|Bat-R_sinicus|Yunnan-Kunming|KT444582|2013-07-21': {'G': 0.20835259161439418, 'C': 0.20066028392208649, 'A': 0.28451634202707166, 'T': 0.30647078243644765}, 'Rs4874|Bat-R_sinicus|Yunnan-Kunming|KY417150|2013-07-21': {'G': 0.20830721520240178, 'C': 0.2005542542311372, 'A': 0.2850450331562799, 'T': 0.3060934974101811}, 'YN2018B|Bat-R_affinis|Yunnan|MK211376|2016-09': {'G': 0.20759518773135907, 'C': 0.20081967213114754, 'A': 0.2848360655737705, 'T': 0.3067490745637229}, 'Rs7327|Bat-R_sinicus|Yunnan--Kunming|KY417151|2014-10-24': {'G': 0.20774078595703963, 'C': 0.20081169366812948, 'A': 0.28528062823770084, 'T': 0.30616689213713005}, 'Rs9401|Bat-R_sinicus|Yunnan-Kunming|KY417152|2015-10-16': {'G': 0.2080352044072693, 'C': 0.20017467835667976, 'A': 0.28506164130471295, 'T': 0.30672847593133795}, 'Rs4084|Bat-R_sinicus|Yunnan-Kunming|KY417144|2012-09-18': {'G': 0.20806180718844475, 'C': 0.20137722539469263, 'A': 0.2850856567013772, 'T': 0.3054753107154854}, 'RsSHC014|Bat-R_sinicus|Yunnan-Kunming|KC881005|2011-04-17': {'G': 0.20807734917917212, 'C': 0.20096015040118173, 'A': 0.2854601000436432, 'T': 0.3054016852989559}, 'Rs3367|Bat-R_sinicus|Yunnan-Kunming|KC881006|2012-03-19': {'G': 0.20784103114930183, 'C': 0.20012083780880774, 'A': 0.28594924812030076, 'T': 0.3060888829215897}, 'WIV1|Bat-R_sinicus|Yunnan-Kunming|KF367457|2012-09': {'G': 0.20766109076511927, 'C': 0.20000659870005608, 'A': 0.2855257514269689, 'T': 0.30680655910785576}, 'YN2018C|Bat-R_affinis|Yunnan-Kunming|MK211377|2016-09': {'G': 0.2088315537741251, 'C': 0.20206136953080264, 'A': 0.28512243591902725, 'T': 0.303984640776045}, 'As6526|Bat-Aselliscus_stoliczkanus|Yunnan-Kunming|KY417142|2014-05-12': {'G': 0.20904962153069806, 'C': 0.2016484440706476, 'A': 0.28518082422203533, 'T': 0.304121110176619}, 'YN2018D|Bat-R_affinis|Yunnan|MK211378|2016-09': {'G': 0.20855260980372686, 'C': 0.20173435276205606, 'A': 0.28540694403071526, 'T': 0.3043060934035018}, 'Rs4081|Bat-R_sinicus|Yunnan-Kunming|KY417143|2012-09-18': {'G': 0.20920614639722943, 'C': 0.20143909081739014, 'A': 0.2854645102720151, 'T': 0.3038902525133654}, 'Rs4255|Bat-R_sinicus|Yunnan-Kunming|KY417149|2013-04-17': {'G': 0.20899035067074606, 'C': 0.20092122516222305, 'A': 0.2856806643580002, 'T': 0.3044077598090307}, 'Rs4237|Bat-R_sinicus|Yunnan-Kunming|KY417147|2013-04-17': {'G': 0.20920614639722943, 'C': 0.2003631350660704, 'A': 0.28583437006153123, 'T': 0.304596348475169}, 'Rs4247|Bat-R_sinicus|Yunnan-Kunming|KY417148|2013-04-17': {'G': 0.20902397202703157, 'C': 0.200786739737081, 'A': 0.2857142857142857, 'T': 0.3044750025216017}, 'Rs672|Bat-R_sinicus|Guizhou|FJ588686|2006-09': {'G': 0.20940156233869026, 'C': 0.20152104339447333, 'A': 0.2844557624143983, 'T': 0.3046216318524381}, 'YN2018A|Bat-R_affinis|Yunnan|MK211375|2016-09': {'G': 0.20883561182571217, 'C': 0.20166341167755406, 'A': 0.2850360293622466, 'T': 0.3044649471344872}, 'YN2013|Bat-R_sinicus|Yunnan|KJ473816|2010-12': {'G': 0.21045226820396679, 'C': 0.20056962459680186, 'A': 0.28241026696863636, 'T': 0.306567840230595}, 'Anlong-103|Bat-R_sinicus|Guizhou-Anlong|KY770858|2013': {'G': 0.20954594448935596, 'C': 0.20112503368364323, 'A': 0.283447857720291, 'T': 0.3058811641067098}, 'Anlong-112|Bat-R_sinicus|Guizhou-Anlong|KY770859|2013': {'G': 0.20969427309805508, 'C': 0.20109886405770722, 'A': 0.2830754710621229, 'T': 0.3061313917821148}, 'HSZ-Cc|SARS-CoV-1|Guangzhou|AY394995|2002': {'G': 0.20792877540735763, 'C': 0.1998656139761465, 'A': 0.28509994960524104, 'T': 0.3071056610112548}, 'YNLF_31C|Bat-R_Ferrumequinum|Yunnan-Lufeng|KP886808|2013-05-23': {'G': 0.20828987652659556, 'C': 0.19910507014769707, 'A': 0.2850654375399522, 'T': 0.3075396157857551}, 'YNLF_34C|Bat-R_Ferrumequinum|Yunnan-Lufeng|KP886809|2013-05-23': {'G': 0.20832352050600544, 'C': 0.19910507014769707, 'A': 0.285132725498772, 'T': 0.3074386838475255}, 'F46|Bat-R_pusillus|Yunnan|KU973692|2012': {'G': 0.21031559114460668, 'C': 0.20136599152143195, 'A': 0.2837965143664625, 'T': 0.3045219029674988}, 'SC2018|Bat-R_spp|Sichuan|MK211374|2016-10': {'G': 0.2096937398812736, 'C': 0.19913653534808418, 'A': 0.28332433890987585, 'T': 0.3078453858607663}, 'LYRa11|Bat-R_affinis|Yunnan-Baoshan|KF569996|2011': {'G': 0.20848850863949, 'C': 0.19832242912263043, 'A': 0.2843482637141419, 'T': 0.30884079852373764}, 'Yunnan2011|Bat-Chaerephon_plicata|Yunnan|JX993988|2011-11': {'G': 0.2092557381502105, 'C': 0.19947711530626103, 'A': 0.2844628548146136, 'T': 0.3068042917289148}, 'Longquan_140|Bat-R_monoceros|China|KF294457|2012': {'G': 0.20996765062676911, 'C': 0.1983420946219167, 'A': 0.28531473244372557, 'T': 0.30637552230758863}, 'HKU3-1|Bat-R_sinicus|Hong_Kong|DQ022305|2005-02-17': {'G': 0.21131593110871905, 'C': 0.19987890204520992, 'A': 0.2841428955866523, 'T': 0.30466227125941875}, 'HKU3-3|Bat-R_sinicus|Hong_Kong|DQ084200|2005-03-17': {'G': 0.21133586887011543, 'C': 0.19989229578270673, 'A': 0.28420450338258557, 'T': 0.30456733196459224}, 'HKU3-2|Bat-R_sinicus|Hong_Kong|DQ084199|2005-02-24': {'G': 0.2116751440024253, 'C': 0.19995284131101154, 'A': 0.28342372082056116, 'T': 0.304948293866002}, 'HKU3-4|Bat-R_sinicus|Hong_Kong|GQ153539|2005-07-20': {'G': 0.21172232695933207, 'C': 0.19997306760032318, 'A': 0.2832951791004579, 'T': 0.30500942633988687}, 'HKU3-5|Bat-R_sinicus|Hong_Kong|GQ153540|2005-09-20': {'G': 0.211553999461352, 'C': 0.20010772959870723, 'A': 0.2833961755992459, 'T': 0.30494209534069483}, 'HKU3-6|Bat-R_sinicus|Hong_Kong|GQ153541|2005-12-16': {'G': 0.21168866145973605, 'C': 0.20007406409911124, 'A': 0.28332884460005386, 'T': 0.30490842984109884}, 'HKU3-10|Bat-R_sinicus|Hong_Kong|GQ153545|2006-10-28': {'G': 0.2115844418252231, 'C': 0.1999663242970197, 'A': 0.2832800134702812, 'T': 0.305169220407476}, 'HKU3-9|Bat-R_sinicus|Hong_Kong|GQ153544|2006-10-28': {'G': 0.2115844418252231, 'C': 0.1999326485940394, 'A': 0.2832800134702812, 'T': 0.3052028961104563}, 'HKU3-11|Bat-R_sinicus|Hong_King|GQ153546|2007-03-07': {'G': 0.2115844418252231, 'C': 0.2, 'A': 0.2832800134702812, 'T': 0.3051355447044957}, 'HKU3-13|Bat-R_sinicus|Hong_Kong|GQ153548|2007-11-15': {'G': 0.21147690130404018, 'C': 0.2001886983185632, 'A': 0.28328335074299965, 'T': 0.305051049634397}, 'HKU3-12|Bat-R_sinicus|Hong_Kong|GQ153547|2007-05-15': {'G': 0.21098168596821976, 'C': 0.19980474010234311, 'A': 0.28373283059520604, 'T': 0.3054807433342311}, 'HKU3-7|Bat-R_sinicus|Guangdong|GQ153542|2006-02-15': {'G': 0.21187239197738592, 'C': 0.1991519720016153, 'A': 0.2832144299367344, 'T': 0.3057612060842644}, 'HKU3-8|Bat-R_sinicus|Guangdong|GQ153543|2006-02-15': {'G': 0.21192008355513628, 'C': 0.1989151308918163, 'A': 0.2831104073312894, 'T': 0.306054378221758}, 'CoVZC45|Bat-R_sinicus|Zhoushan-Dinghai|MG772933|2017-02': {'G': 0.20199986578082008, 'C': 0.1870344272196497, 'A': 0.2932689081269713, 'T': 0.31769679887255886}, 'CoVZXC21|Bat-R_sinicus|Zhoushan-Dinghai|MG772934|2015-07': {'G': 0.20099556033902866, 'C': 0.18723933808690973, 'A': 0.2937575676039284, 'T': 0.3180075339701332}, 'Wuhan-Hu-1|SARS-CoV-2|Wuhan|MN908947|2019-12': {'G': 0.19606728421897468, 'C': 0.18366050229074005, 'A': 0.29943483931378123, 'T': 0.32083737417650404}, 'BtKY72|Bat-R_spp|Kenya|KY352407|2007-10': {'G': 0.2067363530778165, 'C': 0.18542050966728155, 'A': 0.2853726856596297, 'T': 0.32247045159527227}, 'BM48-31|Bat-R_blasii|Bulgaria|NC_014470|2008-04': {'G': 0.2101038393223118, 'C': 0.19439131028829076, 'A': 0.27770187184041534, 'T': 0.3178029785489821}, 'RaTG13|Bat-R_affinis|Yunnan|EPI_ISL_402131|2013-07-24': {'G': 0.1958130966337297, 'C': 0.18455870038519512, 'A': 0.29904538603249037, 'T': 0.3205828169485848}, 'P4L|pangolin|Guangxi|EPI_ISL_410538|2017': {'G': 0.19714151513118164, 'C': 0.18804938602965846, 'A': 0.29849694692343826, 'T': 0.31627860162383414}, 'P5L|pangolin|Guangxi||EPI_ISL_410540|2017': {'G': 0.19717506542306917, 'C': 0.18788163457022075, 'A': 0.29859759779910083, 'T': 0.3163457022076092}, 'P5E|pangolin|Guangxi|EPI_ISL_410541|2017': {'G': 0.1972757162987318, 'C': 0.18801583573777092, 'A': 0.2985640475072133, 'T': 0.3160101992887338}, 'P1E|pangolin|Guangxi|EPI_ISL_410539|2017': {'G': 0.19717459145666252, 'C': 0.1878796013556592, 'A': 0.2987148082279118, 'T': 0.3162309989597665}, 'P2V|pangolin|Guangxi|EPI_ISL_410542|2017': {'G': 0.19714717234435308, 'C': 0.1878503104547743, 'A': 0.2985064608155731, 'T': 0.31616042960228224}, 'Pangolin-CoV|pangolin|Guandong|EPI_ISL_410721|2020-02-16': {'G': 0.1965898369052433, 'C': 0.18550343712090578, 'A': 0.3002763175630139, 'T': 0.3172934357730152}, 'RmYN02|EPI_ISL_412977|China|Yunnan|Xishuangbanna|2019-06-25': {'G': 0.19760035051059957, 'C': 0.18479323244919282, 'A': 0.29746216844730544, 'T': 0.32014424859290214}}


####    CUSTOM DICTIONARIES    ####

#prints out tab-delimited header
print('sample\tacc\tpos1\tpos2')

#for each virus 
for theid in list(full_ntcont):

    #take 10 samples per virus 
    for j in range(10):

        #define the dictionary with single nucleotide estimates 
        ntcont = full_ntcont[theid]

        all_nt = []

        #for each nucleotide
        for n in list(ntcont):
            #determine number of occurences of this nucleotide given its 
            #whole-genome frequency and total number of nucleotides used here
            i = int(ntcont[n] * ntsamples) 
            #add that number of this nucleotide in the list
            for ii in range(i):
                all_nt.append(n)

        #shuffle the nucleotides randomly in the list
        all_nt_sh = random.sample(all_nt, len(all_nt))

        all_cods = []

        #make sure number of nucleotides in the list is a multiple of 3
        remainder = len(all_nt_sh) % 3
        if remainder != 0:
            all_nt_sh = all_nt_sh[:-remainder]
        #put the randomised nucleotides into codons (groups of 3)
        for i in range(0, len(all_nt_sh), 3):
            all_cods.append(str(all_nt_sh[i] + all_nt_sh[i+1] + all_nt_sh[i+2])) 


        syco_new = {}

        #populate an amino acid dictionary : synonymous codons dictionary
        #with the generated codons
        for a in list(syco_orig):
            syco_new.update({a:[]})
            for c in syco_orig[a]:
                for ac in all_cods:
                    if ac == c:
                        syco_new[a].append(ac)
                
                
                

        CpGpos1 = {}


        CpGpos2 = {}


        #for each amino acid
        for a in syco_new:
        
            #store the position 1 and 2 dinucleotides for each generated codon

            dionelist = []
            
            ditwolist = []
            
            for cod in syco_new[a]:
            
                dione = cod[:-1]
                
                ditwo = cod[1:]
                
                dionelist.append(dione)
                
                ditwolist.append(ditwo)


            #calculate the proportion of CpG dinucleotides in the generated
            #synonymous codons for each amino acid
            
            if 'CG' in dionelist:
                
                numCpG1 = dionelist.count('CG')
                
                prob1 = numCpG1/len(dionelist)
                
                CpGpos1.update({a:prob1})
                
            if 'CG' in ditwolist:
                
                numCpG2 = ditwolist.count('CG')
                
                prob2 = numCpG2/len(ditwolist)
                
                CpGpos2.update({a:prob2})  


        #print this sample's tab-delimited row
        print(str(j+1) + '\t' + theid + '\t' + str(CpGpos1) + '\t' + str(CpGpos2))
