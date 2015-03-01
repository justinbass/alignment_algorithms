################################################################################
# Monotonic Alignment Sequential-Match/Mismatch-Gain Algorithm (MASMGA)
################################################################################

import math
import random
from functools import partial

letters_dna = ['A','T','G','C']
E = 2.7182818

################################################################################
# Functions and Random Distributions
################################################################################

#If c is 0, then str_len = 1 will return 0, which may be undesirable, as this
#   implies that the first indel/match will be ignored. I suggest a value of c=1
def logarithmic(str_len,a,b,c):
    if str_len + c <= 0:
        return 0
    return a + b*math.log(str_len + c, E)

def linear(str_len):
    return str_len

def affinelog(str_len,a,b,c,d):
    if str_len + d <= 0:
        return 0
    return a + b*str_len + c*math.log(str_len + d, E)

def affine(str_len,a,b):
    return a + b*str_len

def subsubquadratic(str_len,a,b):
    return a + b*math.pow(str_len, 1.1)

def subquadratic(str_len,a,b):
    return a + b*math.pow(str_len, 1.5)

def quadratic(str_len,a,b):
    return a + b*math.pow(str_len, 2)

def cubic(str_len,a,b):
    return a + b*math.pow(str_len, 3)

def exponential(str_len,a,b):
    return a + b*math.pow(E, str_len)

#Empirically measured indel size distribution
def empirical_indel_size_dist():
    return 1.819*math.pow(random.uniform(0.00001,1.061),-0.654)

#Get a bool with uniform frequency
def uniform_bool(freq):
    if freq >= 1.0:
        return True

    if freq <= 0.0:
        return False

    return random.uniform(0,1.0/freq) < 1.0

################################################################################
# Sequence Creation and Modification Functions
################################################################################

#Get a random sequence of letters
def random_seq(letters,length):
    ret_str = ""
    for x in range(0,max(length,1)):
        ret_str += letters[random.randint(0,len(letters)-1)]
    return ret_str

#Introduce SNPs into a string in proportion to accuracy
def add_snps(instr, freq, letters):
    num = 0
    for i in range(0,len(instr)):
        if uniform_bool(freq):
            ins_list = list(instr)
            new_let = letters[random.randint(0,len(letters)-1)]

            #Only increase SNP count if the letter has changed
            if new_let != ins_list[i]:
                num += 1

            ins_list[i] = new_let
            instr = "".join(ins_list)

    return (instr, num)

#Introduce random insertions into a string with a frequency of an insertion
#   starting, and a distribution of the length of insertion once started
def add_random_insertion(instr, freq, dist_length):
    num = 0
    ilen = 0
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(dist_length()))
            instrlist.insert(i,random_seq(letters_dna,length))
            num += 1
            ilen += length

    return ("".join(instrlist),num,ilen)

#Introduce random deletions into a string with a frequency of an deletion
#   starting, and a distribution of the length of deletion once started
def add_random_deletion(instr, freq, dist_length):
    num = 0
    ilen = 0
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(dist_length()))

            #j is only used for counting up. We want to remove at position i to
            #   avoid deleting off the end of the list.
            for j in range(0, min(length,len(instrlist)-i)):
                instrlist.pop(i)

            num += 1
            ilen += min(length,len(instrlist)-i)

    return ("".join(instrlist),num,ilen)

#Introduce random Copy-Number Variations
#   freq: Freqeuncy of a CNV starting
#   dist_length: Distribution of length of copied string
#   dist_copies: Distribution of number of copies
def add_random_cnv(instr, freq, dist_length, dist_copies):
    num = 0
    clen = 0
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(dist_length()))

            #Get string to be copied
            copy_string = ""
            for j in range(0, min(length,len(instrlist)-i)):
                copy_string += instrlist[i+j]
            copy_string = list(copy_string)

            #Reverse so the list is inserted backwards, e.g. reading forwards
            #TODO: If this line is removed, copies will be backwards. This might
            #   be biologically plausible with some small frequency.
            copy_string.reverse()

            #Insert the copy_string a number of times according to dist_copies()
            copies = dist_copies()
            for j in range(0, copies):
                for copy_letter in copy_string:
                    instrlist.insert(i,copy_letter)

            num += 1
            clen += length*copies

    return ("".join(instrlist),num,clen)

################################################################################
# The Main MASMGA Function
################################################################################

def monotonicAlign(seqA,seqB,match_fun,mismatch_fun,gappen_fun):
    #Get score matrix
    scoremat = list()

    #Three lists for keeping track of continuous matches and indels
    Dmat = list()
    Dmis = list()
    Vind = list()
    Hind = list()

    for i in range(len(seqB)+1):
        scoremat.append([0.0]*(len(seqA)+1))
        Dmat.append([0.0]*(len(seqA)+1))
        Dmis.append([0.0]*(len(seqA)+1))
        Vind.append([0.0]*(len(seqA)+1))
        Hind.append([0.0]*(len(seqA)+1))

    for column in range(1,len(seqA)+1):
        Hind[0][column] = Hind[0][column-1] + 1.0
        Hind_score = gappen_fun(Hind[0][column]) - gappen_fun(Hind[0][column-1])
        scoremat[0][column] = scoremat[0][column-1] - Hind_score
        
    for line in range(1,len(seqB)+1):
        Vind[line][0] = Vind[line-1][0] + 1.0
        Vind_score = gappen_fun(Vind[line][0]) - gappen_fun(Vind[line-1][0])
        scoremat[line][0] = scoremat[line-1][0] - Vind_score

    for column in range(1,len(seqA)+1):
        for line in range(1,len(seqB)+1):

            if seqA[column-1] == seqB[line-1]:
                Dmis[line][column] = 0.0

                length = Dmat[line-1][column-1] + 1.0
                Dmat[line][column] = length
                match_score = match_fun(length) - match_fun(length - 1)
            else:
                Dmat[line][column] = 0.0

                length = Dmis[line-1][column-1] + 1.0
                Dmis[line][column] = length
                match_score = -(mismatch_fun(length) - mismatch_fun(length - 1))

            #Set the Vind and Hind matrices to increment for now
            #   if match/mismatch is more optimal, these will be set to 0 later
            Hind[line][column] = Hind[line][column-1] + 1.0
            Vind[line][column] = Vind[line-1][column] + 1.0

            Hind_score = gappen_fun(Hind[line][column]) - gappen_fun(Hind[line][column-1])
            Vind_score = gappen_fun(Vind[line][column]) - gappen_fun(Vind[line-1][column])

            diagonal = scoremat[line-1][column-1] + match_score
            left = scoremat[line][column-1] - Hind_score
            up = scoremat[line-1][column] - Vind_score

            maximum = [diagonal, left, up]
            m = max(maximum)
            scoremat[line][column] = m

            #up == left == diagonal OR diagonal (lazy eval)
            if (len(set(maximum))==1) or (m == diagonal):
                Vind[line][column] = 0.0
                Hind[line][column] = 0.0

            #left
            elif m == left:
                Dmat[line][column] = 0.0
                Dmis[line][column] = 0.0

            #up
            else:
                Dmat[line][column] = 0.0
                Dmis[line][column] = 0.0
            

    if VERBOSE:
        print 'Score:'
        for m in scoremat:
            for i in range(0,len(m)):
                m[i] = round(m[i],2)
            print m
        print ''

        print 'Match Length:'
        for m in Dmat:
            print m
        print ''

        print 'Mismatch Length:'
        for m in Dmis:
            print m
        print ''

        print 'Horizontal Indel Length:'
        for m in Hind:
            print m
        print ''

        print 'Vertical Indel Length:'
        for m in Vind:
            print m
        print ''

    #Traceback:
    a = len(seqA) #column
    b = len(seqB) #line

    alignmentA=""
    alignmentB=""

    if VERBOSE:
        print ''
        print 'Traceback:'
        print 's1 s2'

    while (a != 0) or (b != 0):

        if VERBOSE:
            print a,b

        if b > 0:
            up = scoremat[b-1][a]
        else:
            up = -float('inf')

        if a > 0:
            left = scoremat[b][a-1]
        else:
            left = -float('inf')

        if a > 0 and b > 0:
            diagonal = scoremat[b-1][a-1]
        else:
            diagonal = -float('inf')

        maximum = [diagonal,left,up]
        m = max(maximum)

        #up == left == diagonal OR diagonal (lazy eval)
        if (len(set(maximum))==1) or (m == diagonal):
            a-=1
            b-=1
            alignmentA=seqA[a]+alignmentA
            alignmentB=seqB[b]+alignmentB

        #left
        elif m == left:
            a-=1
            alignmentA=seqA[a]+alignmentA
            alignmentB="-"+alignmentB

        #up
        else:
            b-=1
            alignmentA="-"+alignmentA
            alignmentB=seqB[b]+alignmentB

    #score = scoremat[-1][-1]

    return [alignmentA,alignmentB]

################################################################################
# Testing and Data Collection
################################################################################

#Return the alignment score, matches/(matches+mismatches+indels)
def score_alignment(align1,align2,letters):
    if len(align1) != len(align2):
        return -1

    matches = 0.0
    mismatches = 0.0
    indels = 0.0
    for i in range(0,len(align1)):

        if align1[i] == align2[i] and align1[i] in letters and align2[i] in letters:
            matches += 1.0
        elif align1[i] != align2[i] and align1[i] in letters and align2[i] in letters:
            mismatches += 1.0
        else:
            indels += 1.0

    return [matches, mismatches, indels]

#Get the indel distribution as a dictionary: key=length, value=frequency
def get_indel_dist(align1,align2):
    indel_dist = dict()

    for align in [align1,align2]:
        #Set a fake last letter to allow parsing of the real last letter
        align += ' '

        indel_len = 0
        for c in align:
            if c == '-':
                indel_len += 1
            else:
                if indel_len > 0:
                    if not indel_len in indel_dist:
                        indel_dist[indel_len] = 0

                    indel_dist[indel_len] += 1

                    indel_len = 0

    return indel_dist

VERBOSE = False
PRINT_IN_OUT = False

#funs = [logarithmic,linear,affinelog,affine,subsubquadratic,subquadratic,
#        quadratic,cubic,exponential]

gap_funs = list()

for aa,bb in [(aa,bb) for aa in range(0,5) for bb in range(0,5)]:
    gap_funs.append(partial(affine,a=aa,b=bb))
'''
for aa,bb,cc in [(aa,bb,cc) for aa in range(0,5) for bb in range(0,5) for cc in range(0,5)]:
    gap_funs.append(partial(logarithmic,a=aa,b=bb,c=cc))
for aa,bb,cc,dd in [(aa,bb,cc,dd) for aa in range(0,5) for bb in range(0,5) for cc in range(0,5) for dd in range(0,5)]:
    gap_funs.append(partial(affinelog,a=aa,b=bb,c=cc,d=dd))
for aa,bb in [(aa,bb) for aa in range(0,5) for bb in range(0,5)]:
    gap_funs.append(partial(subsubquadratic,a=aa,b=bb))
for aa,bb in [(aa,bb) for aa in range(0,5) for bb in range(0,5)]:
    gap_funs.append(partial(subquadratic,a=aa,b=bb))
for aa,bb in [(aa,bb) for aa in range(0,5) for bb in range(0,5)]:
    gap_funs.append(partial(quadratic,a=aa,b=bb))
'''

mat_funs = [partial(linear)]
mismat_funs = [partial(linear)]

trials_grid = [(mat_fun, mismat_fun, gap_fun) for mat_fun in mat_funs
                                              for mismat_fun in mismat_funs
                                              for gap_fun in gap_funs]

datacsv = open('data.csv','w')

write_head = 'trials,len1,len2,\
actual_matches,actual_mismatches,actual_indels,actual_indel_mean,\
total_matches,total_mismatches,total_indels,total_indelmean,\
ma_name,ma,mb,mc,md,mma_name,mma,mmb,mmc,mmd,gf_name,ga,gb,gc,gd'

print write_head
datacsv.write(write_head + '\n')

for mat_fun, mismat_fun, gap_fun in trials_grid:
    trials = 100
    str_length = 100
    total_matches = 0
    total_mismatches = 0
    total_indels = 0
    total_indelmean = 0
    len1 = 0
    len2 = 0
    tinsnum = 0
    tinslen = 0
    tdelnum = 0
    tdellen = 0
    tcnvnum = 0
    tcnvlen = 0
    tsnpnum = 0
    
    for t in range(trials):
        str1 = random_seq(letters_dna, str_length)
        str2 = str1

        (str2,insnum,inslen) = add_random_insertion(str2, 0.02, empirical_indel_size_dist)
        (str2,delnum,dellen) = add_random_deletion(str2, 0.02, empirical_indel_size_dist)
        (str2,cnvnum,cnvlen) = add_random_cnv(str2, 0.02, partial(linear,str_len=4),partial(linear,str_len=1))
        (str2,snpnum) = add_snps(str2, 0.02, letters_dna)

        [align1,align2] = monotonicAlign(str1,str2,
                                         mat_fun,
                                         mismat_fun,
                                         gap_fun)

        [matches, mismatches, indels] = score_alignment(align1,align2,letters_dna)
        total_matches += matches
        total_mismatches += mismatches
        total_indels += indels
        len1 += len(str1)
        len2 += len(str2)
        tinsnum += insnum
        tinslen += inslen
        tdelnum += delnum
        tdellen += dellen
        tcnvnum += cnvnum
        tcnvlen += cnvlen
        tsnpnum += snpnum

        def get_indel_mean(idist):
            indel_mean = 0.0
            indel_num = 0.0
            for i in idist:
                indel_mean += i*idist[i]
                indel_num += idist[i]

            if indel_num != 0:
                indel_mean = indel_mean/indel_num
            else:
                indel_mean = 0

            return indel_mean

        total_indelmean += get_indel_mean(get_indel_dist(align1,align2))

        if PRINT_IN_OUT:
            print 'Input seqs:\n'+str1+'\n'+str2+'\n'
            print 'Aligned seqs (score '+str(score)+'):\n'+align1+'\n'+align2
            print 'Indel Distribution: ' + str(indel_dist) + '\n'

    ma_name = mat_fun.func.__name__
    args = mat_fun.keywords
    ma = args['a'] if args and 'a' in args else ''
    mb = args['b'] if args and 'b' in args else ''
    mc = args['c'] if args and 'c' in args else ''
    md = args['d'] if args and 'd' in args else ''

    mma_name = mismat_fun.func.__name__
    args = mismat_fun.keywords
    mma = args['a'] if args and 'a' in args else ''
    mmb = args['b'] if args and 'b' in args else ''
    mmc = args['c'] if args and 'c' in args else ''
    mmd = args['d'] if args and 'd' in args else ''

    gf_name = gap_fun.func.__name__
    args = gap_fun.keywords
    ga = args['a'] if args and 'a' in args else ''
    gb = args['b'] if args and 'b' in args else ''
    gc = args['c'] if args and 'c' in args else ''
    gd = args['d'] if args and 'd' in args else ''

    total_matches /= trials
    total_mismatches /= trials
    total_indels /= trials
    len1 /= trials
    len2 /= trials
    tinsnum /= trials
    tinslen /= trials
    tdelnum /= trials
    tdellen /= trials
    tcnvnum /= trials
    tcnvlen /= trials
    tsnpnum /= trials
    total_indelmean /= trials

    actual_matches = len1 - tdellen
    actual_mismatches = tsnpnum
    actual_indels = tinsnum + tdelnum + tcnvnum
    actual_indel_mean = (tinslen+tdellen+tcnvlen)/actual_indels

    write_vars = [trials,len1,len2,
    actual_matches,actual_mismatches,actual_indels,actual_indel_mean,
    total_matches,total_mismatches,total_indels,total_indelmean,
    ma_name,ma,mb,mc,md, mma_name,mma,mmb,mmc,mmd, gf_name,ga,gb,gc,gd]
    write_str = str(write_vars[0])
    for i in range(1,len(write_vars)):
        write_str += ',' + str(write_vars[i])

    print write_str
    datacsv.write(write_str + '\n')
    datacsv.flush()

