################################################################################
# Non-Gap Match Gain and Mismatch Penalty Alignment Algorithm (NGMA, or ENIGMA)
################################################################################

import math
import random
import time,datetime
from functools import partial

VERBOSE = False
PRINT_IN_OUT = False
PRINT_ETA = True

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

        #Do first column/line separate to account for constant shifts in the gappen
        if(Hind[0][column] == 1.0):
            Hind_score = gappen_fun(Hind[0][column])
        else:
            Hind_score = gappen_fun(Hind[0][column]) - gappen_fun(Hind[0][column-1])

        scoremat[0][column] = scoremat[0][column-1] - Hind_score
        
    for line in range(1,len(seqB)+1):
        Vind[line][0] = Vind[line-1][0] + 1.0

        if(Vind[line][0] == 1.0):
            Vind_score = gappen_fun(Vind[line][0])
        else:
            Vind_score = gappen_fun(Vind[line][0]) - gappen_fun(Vind[line-1][0])

        scoremat[line][0] = scoremat[line-1][0] - Vind_score

    for column in range(1,len(seqA)+1):
        for line in range(1,len(seqB)+1):

            if seqA[column-1] == seqB[line-1]:
                Dmis[line][column] = 0.0

                length = Dmat[line-1][column-1] + 1.0
                Dmat[line][column] = length
                if length == 1.0:
                    match_score = match_fun(1.0)
                else:
                    match_score = match_fun(length) - match_fun(length - 1)
            else:
                Dmat[line][column] = 0.0

                length = Dmis[line-1][column-1] + 1.0
                Dmis[line][column] = length
                if length == 1.0:
                    match_score = -(mismatch_fun(1.0))
                else:
                    match_score = -(mismatch_fun(length) - mismatch_fun(length - 1))

            #Set the Vind and Hind matrices to increment for now
            #   if match/mismatch is more optimal, these will be set to 0 later
            Hind[line][column] = Hind[line][column-1] + 1.0
            Vind[line][column] = Vind[line-1][column] + 1.0

            if Hind[line][column] == 1.0:
                Hind_score = gappen_fun(1.0)
            else:
                Hind_score = gappen_fun(Hind[line][column]) - gappen_fun(Hind[line][column-1])

            if Vind[line][column] == 1.0:
                Vind_score = gappen_fun(1.0)
            else:
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

gap_funs = list()

maxi = 2

gap_funs.append(partial(affine,a=0,b=0))
for aa,bb in [(aa,bb) for aa in range(0,maxi) for bb in range(1,maxi)]:
    gap_funs.append(partial(affine,a=aa,b=bb))
for aa,bb,cc in [(aa,bb,cc) for aa in range(0,maxi) for bb in range(1,maxi) for cc in range(1,maxi)]:
    gap_funs.append(partial(logarithmic,a=aa,b=bb,c=cc))
for aa,bb,cc,dd in [(aa,bb,cc,dd) for aa in range(0,maxi) for bb in range(1,maxi) for cc in range(1,maxi) for dd in range(1,maxi)]:
    gap_funs.append(partial(affinelog,a=aa,b=bb,c=cc,d=dd))
for aa,bb in [(aa,bb) for aa in range(0,maxi) for bb in range(1,maxi)]:
    gap_funs.append(partial(subquadratic,a=aa,b=bb))
for aa,bb in [(aa,bb) for aa in range(0,maxi) for bb in range(1,maxi)]:
    gap_funs.append(partial(quadratic,a=aa,b=bb))

mat_funs = [partial(linear)]
mismat_funs = [partial(linear)]

trials_grid = [(mat_fun, mismat_fun, gap_fun) for mat_fun in gap_funs
                                              for mismat_fun in gap_funs
                                              for gap_fun in gap_funs]

datacsv = open('data.csv','w')

write_head = 'trials,len1,len2,\
actual_matches,actual_mismatches,actual_indels,actual_indel_mean,\
meas_matches,meas_mismatches,meas_indels,meas_indelmean,\
ma_name,ma,mb,mc,md,mma_name,mma,mmb,mmc,mmd,gf_name,ga,gb,gc,gd,\
total_score_mean,total_score_std'

print write_head
datacsv.write(write_head + '\n')

trials = 1000
str_length = 100

w = len(trials_grid)
l = trials
len1 = list()
len2 = list()

actual_matches = list()
actual_mismatches = list()
actual_indels = list()
actual_indelmean = list()

meas_matches = list()
meas_mismatches = list()
meas_indels = list()
meas_indelmean = list()
total_score = list()

for t in range(w):
    meas_matches.append(list())
    meas_mismatches.append(list())
    meas_indels.append(list())
    meas_indelmean.append(list())
    total_score.append(list())

if PRINT_ETA:
    time_begin = time.time()
    print 'BEGIN:',datetime.datetime.now().time()

for t in range(trials):
    str1 = random_seq(letters_dna, str_length)
        
    str2 = str1
    (str2,insnum,inslen) = add_random_insertion(str2, 0.02, empirical_indel_size_dist)
    (str2,delnum,dellen) = add_random_deletion(str2, 0.02, empirical_indel_size_dist)
    (str2,cnvnum,cnvlen) = add_random_cnv(str2, 0.02, partial(linear,str_len=4),partial(linear,str_len=1))
    (str2,snpnum) = add_snps(str2, 0.02, letters_dna)

    len1.append(len(str1))
    len2.append(len(str2))

    actual_matches.append(len1[t] - dellen)
    actual_mismatches.append(snpnum)
    actual_indels.append(insnum + delnum + cnvnum)
    if actual_indels[t] != 0:
        actual_indelmean.append((inslen+dellen+cnvlen)/actual_indels[t])
    else:
        actual_indelmean.append(0)

    i = 0
    for mat_fun, mismat_fun, gap_fun in trials_grid:
        [align1,align2] = monotonicAlign(str1,str2,
                                         mat_fun,
                                         mismat_fun,
                                         gap_fun)

        [matches, mismatches, indels] = score_alignment(align1,align2,letters_dna)
        meas_matches[i].append(matches)
        meas_mismatches[i].append(mismatches)
        meas_indels[i].append(indels)
        
        def get_indel_mean(idist):
            indel_mean = 0.0
            indel_num = 0.0
            for d in idist:
                indel_mean += d*idist[d]
                indel_num += idist[d]

            if indel_num != 0:
                indel_mean = indel_mean/indel_num
            else:
                indel_mean = 0

            return indel_mean

        meas_indelmean[i].append(get_indel_mean(get_indel_dist(align1,align2)))

        def diff_perc(a,b):
            if a == 0 or b == 0:
                return 0
            return min(a/b,b/a)

        ''' #This scoring method means that a single 0 will cause the total result to be 0
        match_score = diff_perc(meas_matches[i][t],actual_matches[t])
        mismatch_score = diff_perc(,actual_mismatches[t])
        indel_score = diff_perc(meas_indels[i][t],actual_indels[t])
        indelmean_score = diff_perc(meas_indelmean[i][t],actual_indelmean[t])
        ts = match_score * mismatch_score * math.sqrt(indel_score*indelmean_score)
        '''
                
        meas_score = float(meas_matches[i][t])/(meas_matches[i][t]+meas_mismatches[i][t]+meas_indels[i][t])
        actual_score = float(actual_matches[t])/(actual_matches[t]+actual_mismatches[t]+actual_indels[t])
        ts = diff_perc(meas_score,actual_score)
        total_score[i].append(ts)

        i += 1

        if PRINT_ETA:
            if i%20 == 0:
                time_elapsed = round(time.time()-time_begin,2)
                total_runs = len(trials_grid)*trials
                this_run = t*len(trials_grid) + i
                eta = int(time_elapsed*total_runs/this_run - time_elapsed)
                eta_h = int(math.floor(eta/3600))
                eta_m = int(math.floor(eta/60 % 60))
                eta_s = eta%60

                print 'Trial: '+str(t+1)+'/'+str(trials)+' '+\
                'Func: '+str(i)+'/'+str(len(trials_grid))+' '+\
                'Elapsed: '+str(time_elapsed)+'s '+\
                'ETA:',str(eta_h)+'h',str(eta_m)+'m',str(eta_s)+'s'

#Go through each function and print out the data
def get_mean(l):
    if len(l) == 0:
        return 0.0

    return float(sum(l))/len(l)
    
def get_var(l):
    mean = get_mean(l)
    if mean == 0.0:
        return 0.0

    var = 0.0
    for i in l:
        var += (i - mean)*(i - mean)
    var /= len(l)
    var = math.sqrt(var)

    return var

len1 = get_mean(len1)
len2 = get_mean(len2)

actual_matches = get_mean(actual_matches)
actual_mismatches = get_mean(actual_mismatches)
actual_indels = get_mean(actual_indels)
actual_indelmean = get_mean(actual_indelmean)

i = 0
for mat_fun, mismat_fun, gap_fun in trials_grid:
    meas_matches[i] = get_mean(meas_matches[i])
    meas_mismatches[i] = get_mean(meas_mismatches[i])
    meas_indels[i] = get_mean(meas_indels[i])
    meas_indelmean[i] = get_mean(meas_indelmean[i])
    
    total_score_mean = get_mean(total_score[i])
    total_score_var = get_var(total_score[i])

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

    write_vars = [trials,len1,len2,
    actual_matches,actual_mismatches,actual_indels,actual_indelmean,
    meas_matches[i],meas_mismatches[i],meas_indels[i],meas_indelmean[i],
    ma_name,ma,mb,mc,md, mma_name,mma,mmb,mmc,mmd, gf_name,ga,gb,gc,gd,
    total_score_mean,total_score_var]
    write_str = str(write_vars[0])
    for w in range(1,len(write_vars)):
        write_str += ',' + str(write_vars[w])

    #print write_str
    datacsv.write(write_str + '\n')
    datacsv.flush()

    i += 1

if PRINT_ETA:
    print 'DONE:',datetime.datetime.now().time()

