import math
import random
from functools import partial

#Needleman-Wunsch Contiguous-Match-Promotion

letters_dna = ['A','T','G','C']

#Get a random sequence of letters
def random_seq(letters,length):
    if length < 1:
        length = 1

    ret_str = ""

    for x in range(0,length):
        ret_str += letters[random.randint(0,len(letters)-1)]

    return ret_str

#Get a bool with uniform frequency
def uniform_bool(freq):
    if freq >= 1.0:
        return True

    if freq <= 0.0:
        return False

    return random.uniform(0,1.0/freq) < 1.0

#Introduce SNPs into a string in proportion to accuracy
def add_snps(instr, letters, freq):
    for i in range(0,len(instr)):
        if uniform_bool(freq):
            ins_list = list(instr)
            ins_list[i] = letters[random.randint(0,len(letters)-1)]
            instr = "".join(ins_list)
    return instr

#Introduce random insertions into a string with a frequency of an indel starting,
#   and a distribution of the length of indels: length(frequency)
#Current version has hard-coded distribution measured empirically
def add_random_insertion(instr, freq):
    instrlist = list(instr)

    #Must traverse instrlist backwards so that inserts occur in the correct order
    for i in reversed(range(0,len(instrlist))):
        if uniform_bool(freq):
            length = int(math.floor(1.819*math.pow(random.uniform(0.00001,1.061),-0.654)))
            instrlist.insert(i,random_seq(letters_dna,length))

    return "".join(instrlist)

#def add_random_deletion
#def add_random_cnv

##########################
# Needleman Wunsch Alignment
##########################

def logarithmic(str_len,a,b):
    if str_len == 0:
        return 0
    return a + b*math.log(str_len,2.7182818)

#The default Needleman-Wunsch algorithm uses all linear input functions
def linear(str_len):
    return str_len

def affine(str_len,a,b):
    return a + b*str_len

def subsubquadratic(str_len,a,b):
    return a + b*math.pow(str_len,1.1)

def subquadratic(str_len,a,b):
    return a + b*math.pow(str_len,1.5)

def quadratic(str_len,a,b):
    return a + b*math.pow(str_len,2)

def cubic(str_len,a,b):
    return a + b*math.pow(str_len,3)

def exponential(str_len,a,b):
    return a + b*math.pow(2.7182818,str_len)

def needlemanWunsch(seqA,seqB,match_fun,mismatch_fun,gappen_fun):
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

    while (a != 0) and (b != 0):
        up = scoremat[b-1][a]
        left = scoremat[b][a-1]
        diagonal = scoremat[b-1][a-1]

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

    return [alignmentA,alignmentB,scoremat[-1][-1]]

VERBOSE = False

str1 = random_seq(letters_dna, 20)

str2 = add_random_insertion(str1, 0.05)
str2 = add_snps(str2, letters_dna, 0.1)

str1 = 'ACGACGCTAAGATCGGGCCA'
str2 = 'ACCGTGGGGCACTAAACGGCTAGCAACCCTAATGGTTTTGTACATAATAAGACGGACGATCGTAAGCCATGATCGGGCCA'

#str1 = 'GCATGCU'
#str2 = 'GATTACA'

print 'Input sequences:'
print str1
print str2
print ''

[align1,align2,score] = needlemanWunsch(str1,str2,partial(exponential,a=0,b=1),partial(linear),partial(logarithmic,a=0,b=1))

print 'Aligned sequences with score '+str(score)+':'
print align1
print align2
print ''

[align1,align2,score] = needlemanWunsch(str1,str2,partial(quadratic,a=0,b=1),partial(linear),partial(logarithmic,a=0,b=1))

print 'Aligned sequences with score '+str(score)+':'
print align1
print align2


