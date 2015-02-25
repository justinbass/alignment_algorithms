import math

#Needleman-Wunsch Contiguous-Match-Promotion

##########################
# Needleman Wunsch Alignment
##########################

def needlemanWunsch(seqA,seqB,minimumIdentity):
    '''
    def ambiguityCode(sequence1,sequence2):
        iupac={ 'AC': 'M',
                'AG': 'R',
                'AT': 'W',
                'CG': 'S',
                'CT': 'Y',
                'GT': 'K',
                'AN': 'N',
                'GN': 'N',
                'CN': 'N',
                'TN': 'N',
                'NT': 'N',
                'TY': 'C',
                'CY': 'C'}
        consensus=""
        for i in range(len(sequence1)):
            l=[sequence1[i],sequence2[i]]
            if "-" in l:l.remove("-")
            l=list(set(l))
            l.sort()

            if len(l)>1:
                letter="".join(l)
                if letter not in iupac:
                    consensus+="N"
                else:
                    consensus+=iupac[letter]
            else:
                consensus+=l[0]
        return consensus

    def perIdentity(finalA,finalB,per=0):
        for k in range(len(finalA)):
            if finalA[k] == finalB[k]:
                per+=1

        return  per*100.0/len(finalA) if len(finalA)>0 else 0
    '''

    def contig_match_fun(str_len):
        return math.pow(str_len,2)

    #Get initial score matrix
    matrix = list()
    contig_match_length = list()
    for i in range(len(seqB)+1):
        matrix.append([0]*(len(seqA)+1))
        contig_match_length.append([0]*(len(seqA)+1))

    match_score = 1
    mismatch_score = -1
    gap_penalty = -1

    matrix[0][0] = 0.0

    for column in range(1,len(seqA)+1):
        matrix[0][column] = matrix[0][column-1] + gap_penalty
        
    for line in range(1,len(seqB)+1):
        matrix[line][0] = matrix[line-1][0] + gap_penalty

    for column in range(1,len(seqA)+1):
        for line in range(1,len(seqB)+1):
            maximum=[]

            if seqA[column-1] == seqB[line-1]:
                old_len = contig_match_length[line-1][column-1]
                new_len = old_len + match_score
                contig_match_length[line][column] = new_len
                
                match_mismatch_addition = contig_match_fun(new_len) - contig_match_fun(old_len)
            else:
                contig_match_length[line][column] = 0
                match_mismatch_addition = mismatch_score

            maximum.append(matrix[line-1][column-1] + match_mismatch_addition)
            maximum.append(matrix[line][column-1]-1)
            maximum.append(matrix[line-1][column]-1)
            matrix[line][column] = max(maximum)

    for m in matrix:
        print m
    print ''

    for m in contig_match_length:
        print m
    print ''

    #Trace-back, requiring matrix, seqA, and seqB:
    c=len(seqA)
    l=len(seqB)

    alignmentA=""
    alignmentB=""

    while (l!=0) and (c!=0):
        #current=matrix[l][c]
        up=matrix[l-1][c]
        front=matrix[l][c-1]
        diagonal=matrix[l-1][c-1]

        maximo=[diagonal,front,up]

        m=max(maximo)

        #i=maximo.index(m)

        #diagonal
        #up == front == diagonal or i == 0
        if (len(set(maximo))==1) or (m == diagonal):
            l-=1;c-=1
            alignmentA=seqA[c]+alignmentA
            alignmentB=seqB[l]+alignmentB

        #front
        #i == 1
        elif m == front:
            c-=1
            alignmentA=seqA[c]+alignmentA
            alignmentB="-"+alignmentB

        #up
        #i == 2, e.g. m == up
        else:
            l-=1
            alignmentA="-"+alignmentA
            alignmentB=seqB[l]+alignmentB

    '''
    perc=perIdentity(alignmentA,alignmentB)
    consensusSequence=""
    if perc>=minimumIdentity:consensusSequence=ambiguityCode(alignmentA,alignmentB)
    '''

    return [alignmentA,alignmentB]#,perc,consensusSequence]


print needlemanWunsch('GCATGCU','GATTACA',0)

