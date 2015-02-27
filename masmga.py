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

#The default Needleman-Wunsch algorithm uses all linear input functions
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

#TODO: This has not been checked
def empirical_ins_size_dist():
    return 1.414*1.819*math.pow(random.uniform(0.00001,1.061),-0.654)

#TODO: This has not been checked
def empirical_del_size_dist():
    return (1.0/1.414)*1.819*math.pow(random.uniform(0.00001,1.061),-0.654)

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
    for i in range(0,len(instr)):
        if uniform_bool(freq):
            ins_list = list(instr)
            ins_list[i] = letters[random.randint(0,len(letters)-1)]
            instr = "".join(ins_list)
    return instr

#Introduce random insertions into a string with a frequency of an insertion
#   starting, and a distribution of the length of insertion once started
def add_random_insertion(instr, freq, dist_length):
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(dist_length()))
            instrlist.insert(i,random_seq(letters_dna,length))

    return "".join(instrlist)

#Introduce random deletions into a string with a frequency of an deletion
#   starting, and a distribution of the length of deletion once started
def add_random_deletion(instr, freq, dist_length):
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(dist_length()))

            #j is only used for counting up. We want to remove at position i to
            #   avoid deleting off the end of the list.
            for j in range(0, min(length,len(instrlist)-i)):
                instrlist.pop(i)

    return "".join(instrlist)

#Introduce random Copy-Number Variations
#   freq: Freqeuncy of a CNV starting
#   dist_length: Distribution of length of copied string
#   dist_copies: Distribution of number of copies
def add_random_cnv(instr, freq, dist_length, dist_copies):
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

    return "".join(instrlist)

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
    mismatches_indels = 0.0
    for i in range(0,len(align1)):

        if align1[i] == align2[i] and align1[i] in letters and align2[i] in letters:
            matches += 1.0
        else:
            mismatches_indels += 1.0

    return matches / (matches + mismatches_indels)

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

PRINT_IN_OUT = True
VERBOSE = False

str1 = random_seq(letters_dna, 20)
str2 = str1

str2 = add_random_insertion(str2, 0.05, empirical_ins_size_dist)
str2 = add_random_deletion(str2, 0.05, empirical_del_size_dist)
str2 = add_random_cnv(str2, 0.1, partial(linear,str_len=4),partial(linear,str_len=1))
str2 = add_snps(str2, 0.1, letters_dna)

#Test case 1
#str1 = 'GCATGCU'
#str2 = 'GATTACA'

#Test case 2
#str1 = 'ACGACGCTAAGATCGGGCCA'
#str2 = 'ACCGTGGGGCACTAAACGGCTAGCAACCCTAATGGTTTTGTACATAATAAGACGGACGATCGTAAGCCATGATCGGGCCA'

#Test case 3
#str1 = 'AAA'
#str2 = 'ACCGTGGGGCACTAAACGGCTAGCA'

if PRINT_IN_OUT:
    print 'Input sequences:'
    print str1
    print str2
    print ''

[align1,align2] = monotonicAlign(str1,str2,
                                 partial(affine,a=0,b=2),
                                 partial(linear),
                                 partial(affine,a=4,b=2))

score = score_alignment(align1,align2,letters_dna)
indel_dist = get_indel_dist(align1,align2)

if PRINT_IN_OUT:
    print 'Aligned sequences with score '+str(score)+':'
    print align1
    print align2
    print 'Indel Distribution: ' + str(indel_dist)
    print ''


