import random

random.seed()

VERBOSE = False #For printing intermediate steps and debugging
letters = ['A','T','G','C']
ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def get_random_seq(len):
    if len < 1:
        len = 1

    ret_str = ""

    for x in range(0,len):
        ret_str += letters[random.randint(0,3)]

    return ret_str

def split_seq(seq, num, smallest, largest, accuracy):
    if len(seq) < 1:
        return []

    if num < 1:
        num = 1

    if smallest < 1:
        smallest = 1

    if largest < 1:
        largest = 1

    if smallest > largest:
        temp = largest
        largest = smallest
        smallest = temp

    if smallest > len(seq):
        smallest = len(seq)

    if largest > len(seq):
        largest = len(seq)

    read_arr = []

    for n in range(0,num):
        read_size = random.randint(smallest,largest)
        read_start_pos = random.randint(0,len(seq)-read_size)
        read_end_pos = read_start_pos + read_size
        read_arr.append(seq[read_start_pos:read_end_pos])

    return read_arr

def naive_polynomial_align(read_arr):
    #Must be sorted/reversed twice to maintain original order... for some reason
    #TODO
    read_arr = sorted(read_arr,key=len)[::-1]
    read_arr = sorted(read_arr,key=len)[::-1]

    if VERBOSE:
        print 'Length-sorted reads:'
        print str(read_arr)
        print ''

    align_score_arrs = list()
    align_shift_arrs = list()

    #For each pair of reads [loops: N^2. Total: 0(N^2 * L^2)]
    for r1 in range(0,len(read_arr)-1):
        read1 = read_arr[r1]
        align_score_arr = list()
        align_shift_arr = list()

        for r2 in range(r1+1,len(read_arr)):
            read2 = read_arr[r2]

            if read1 == read2:
                align_score_arr.append(len(read1))
                align_shift_arr.append(0)
                continue

            #Find the most probable alignment (largest # of matches) [O(L^2)]
            max_alignment = 0
            corresponding_shift = len(read1) #Nonsense value

            for shift in range(1-len(read2),len(read1)):
                alignment = 0

                overlap1 = read1[max(0,shift):min(len(read1),shift+len(read2))]
                overlap2 = read2[max(0,-shift):min(len(read2),len(read1)-shift)]
                #Check that these 2 lengths are the same?

                for i in range(0,len(overlap1)):
                    if overlap1[i] == overlap2[i]:
                        alignment += 1

                if alignment > max_alignment:
                    max_alignment = alignment
                    corresponding_shift = shift

                #Could add: if max_alignment = min(len(read1),len(read2)): break
                #Since it cannot be higher than the smaller of the lengths

            #Insert max alignment and its degree into the align_arr
            align_score_arr.append(max_alignment)
            align_shift_arr.append(corresponding_shift)

        align_score_arrs.append(align_score_arr)
        align_shift_arrs.append(align_shift_arr)

    #Fill in score and shift arrays, which are triangular
    #   but need to be square [0(N^2)]
    align_score_arrs.append(list())
    align_shift_arrs.append(list())

    for j in range(0,len(align_score_arrs)):
        align_score_arrs[j].insert(0,-1)
        align_shift_arrs[j].insert(0,0)
        for i in range(0,j)[::-1]:
            align_score_arrs[j].insert(0,align_score_arrs[i][j-len(read_arr)])
            align_shift_arrs[j].insert(0,-align_shift_arrs[i][j-len(read_arr)])

    #Now align_shift_arrs[i][j] is the amount needed to shift j to align to i

    if VERBOSE:
        print 'Alignment-score array:'
        for a in align_score_arrs:
            print a
        print ''

    #Get alignments in the order they should be processed based on max score
    align_order = list()

    for rounds in range(0,len(align_score_arrs)*len(align_score_arrs)-1):
        i_max = -1
        j_max = -1
        score_max = -1
        for i in range(0,len(align_score_arrs)):
            for j in range(0,len(align_score_arrs[i])):
                
                if (min(i,j),max(i,j)) in align_order:
                    continue

                if align_score_arrs[i][j] > score_max:
                    score_max = align_score_arrs[i][j]
                    i_max = i
                    j_max = j

        #If the largest score is 0, then it is in a seperate mergelet completely
        #Set it to be outside of the read_arr, as a tag for later
        if score_max == 0:
            first_arr = list()
            second_arr = list()
            for (first,second) in align_order:
                first_arr.append(first)
                second_arr.append(second)

            if not max(i_max,j_max) in first_arr:
                if not max(i_max,j_max) in second_arr:
                    align_order.append((len(read_arr),max(i_max,j_max)))
        else:
            align_order.append((min(i_max,j_max),max(i_max,j_max)))

    '''
        #Remove smaller read from consideration by setting its column/row to -1
        #Using the below methods requires rounds
        #   to be at most len(align_score_arrs)-1
        for i in range(0,len(align_score_arrs[max(i_max,j_max)])):
            align_score_arrs[max(i_max,j_max)][i] = -1
        
        for a in align_score_arrs:
            a[max(i_max,j_max)] = -1
    '''

    if VERBOSE:
        print 'Alignment order:'
        print align_order
        print ''
        print 'Alignment-score\'s corresponding shift array'
        for a in align_shift_arrs:
            print a
        print ''

    #Create a 'mergelet' array, combining/sorting reads based on their shift.
    #   This is necessary in case there are two totally separate contigs, which
    #   then must be returned as a mergelet array, since it can't be a string.
    #   If the returned mergelet array has only 1 mergelet, it can be a
    #   single string. Each mergelet in the array is a separate string/contig.
    #Mergelet array structure: [ [mergelet], [mergelet], ... ]
    #Mergelet structure: [(read1_no, shift1), (read2_no, shift2), ... ] 
    mergelet_arr = list()

    #Insert a read_no into the mergelet given another one that may already exist
    def insert_into_mergelet(mergelet, read_no, shift):
        mergelet_inserted = False

        for j in range(0,len(mergelet)):
            #Each mergelet's reads will be inherently sorted by shift
            if mergelet[j][1] > shift:
                mergelet.insert(j,(read_no,shift))
                mergelet_inserted = True
                break

        #If it hasn't been inserted, it must be bigger than all other shifts
        #   in the mergelet.
        if not mergelet_inserted:
            mergelet.append((read_no,shift))

    if VERBOSE:
        print 'Alignment succession:'

    for (first,second) in align_order:
        #Check through each mergelet to find where to insert (read_no, shift)
        #If the first number doesn't exist, we need to insert a new mergelet
        skip_reads = False
        need_new_mergelet = True
        mergelet_to_insert_into = 0

        for mergelet in mergelet_arr:
            #Create an array of reads already existing in the mergelet_arr
            first_arr = list()
            for read in mergelet:
                first_arr.append(read[0])

            #If first is in array, try to add the second in
            if first in first_arr:
                skip_reads = True

                #If second exists, skip this mergelet
                if second in first_arr:
                    break

                #Find location of first
                first_pos = first_arr.index(first)

                #Get shift from first to second, then add shift of first
                insert_shift = align_shift_arrs[first][second]
                insert_shift += mergelet[first_pos][1]

                insert_into_mergelet(mergelet, second, insert_shift)

                break

            #If second is in array, try to add the first in
            if second in first_arr:
                skip_reads = True

                #If first exists, skip this mergelet
                if first in first_arr:
                    break

                #Find location of second
                second_pos = first_arr.index(second)

                #Get shift from second to first, then add shift of second
                insert_shift = align_shift_arrs[second][first]
                insert_shift += mergelet[second_pos][1]

                insert_into_mergelet(mergelet, first, insert_shift)

                break

        if skip_reads:
            if VERBOSE:
                print str((first,second)) + ': ' + str(mergelet_arr)
            continue

        #A new mergelet must be created
        #If the first is out of the array, this means it was previously
        #   tagged as completely separate from the other reads.
        if first == len(read_arr):
            mergelet_arr.append([(second,0)])
        else:
            mergelet_arr.append([(first,0)])
            second_shift = align_shift_arrs[first][second]

            #Insert before or after first entry
            if second_shift < 0:
                mergelet_arr[-1].insert(0,(second,second_shift))
            else:
                mergelet_arr[-1].append((second,second_shift))

        if VERBOSE:
            print str([first,second]) + ': ' + str(mergelet_arr)

    if VERBOSE:
        print ''
        print 'Merging mergelets:'

    #Merge mergelets if possible
    RESTART = False

    i = 0
    while(True):
        for (read1,shift1) in mergelet_arr[i]:
            for j in range(i+1,len(mergelet_arr)):
                for (read2,shift2) in mergelet_arr[j]:
                    if read1 == read2:
                        while len(mergelet_arr[j]) > 0:
                            insert_read = mergelet_arr[j].pop(0)
                            if (insert_read[0], insert_read[1]+shift1-shift2) not in mergelet_arr[i]:
                                insert_into_mergelet(mergelet_arr[i], insert_read[0], insert_read[1]+shift1-shift2)

                        if VERBOSE:
                            print mergelet_arr

                        RESTART = True

                        break

                if RESTART == True:
                    break

            if RESTART == True:
                break

        #If a match was made, the loop must be completely restart in case values
        #   have been inserted before the iterators in the loops; in this case,
        #   they would never be reached.
        if RESTART == True:
            RESTART = False
            i = 0
            continue

        i += 1
        if i >= len(mergelet_arr):
            break

    while [] in mergelet_arr:
        mergelet_arr.remove([])

    if VERBOSE:
        print '\nFinal mergelet array:\n',mergelet_arr,'\n'

    #Get aligned_array, including aligned reads and predicted final strings
    if VERBOSE:
        print 'Aligned-reads array: '

    aligned_array = list()
    for mergelet in mergelet_arr:
        read_shift_arr = list()
        final_string_arr = list()

        #Find min shift
        min_shift = mergelet[0][1]
        for read in mergelet:
            if read[1] < min_shift:
                min_shift = read[1]

        #Replace read_no with actual read, and normalize shifts to begin at 0
        for j in range(0,len(mergelet)):
            read_shift_arr.append( (read_arr[ mergelet[j][0] ],
                              mergelet[j][1]-min_shift ) )

        #Find last char's shift; maybe unnecessary?
        max_shift = 0
        for read in mergelet:
            if read[1] + len(read) - 1 > min_shift:
                max_shift = read[1] + len(read) - 1

        final_str = ''
        for (read, shift) in read_shift_arr:
            i = 0

            for char in read:
                try:
                    final_string_arr[i+shift]
                except IndexError:
                    final_string_arr.insert(i+shift, list())

                final_string_arr[i+shift].append(char)

                i += 1

        for char_arr in final_string_arr:
            char_histo = dict()
            for char in char_arr:
                if not char in char_histo:
                    char_histo[char] = 0
                else:
                    char_histo[char] += 1

            final_str += max(char_histo.iterkeys(),
                             key=(lambda key: char_histo[key]))

        aligned_array.append((final_str, final_string_arr, read_shift_arr))

        if VERBOSE:
            print str((final_str, final_string_arr, read_shift_arr))
            print ''

    return aligned_array

#Not my own code: from Stack Overflow
def longest_substring(string1, string2):
    answer = ""
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                if (len(match) > len(answer)): answer = match
                match = ""
    return answer

def test_naive_single(rseq,seq_reads):
    PRINT_PASS_AND_DATA = False

    aligned_array = naive_polynomial_align(seq_reads)

    #Substring test
    if PRINT_PASS_AND_DATA:
        print 'data:',seq_reads,'\n'
    for i in range(0,len(aligned_array)):
        if aligned_array[i][0] in rseq:
            if PRINT_PASS_AND_DATA:
                print 'PASSED:', aligned_array[i][0], 'is in', rseq
        else:
            print '*FAILED:', aligned_array[i][0], 'is not in', rseq

    #Common substring test: No two returned strings should have a common
    #   substring (they should have been joined in the algorithm if possible)
    test2_failed = False
    for i in range(0,len(aligned_array)):
        for j in range(i+1,len(aligned_array)):
            substr = longest_substring(aligned_array[i][0],aligned_array[j][0])
            if len(substr) > 0:
                test2_failed = True
                print '*FAILED:', aligned_array[i][0], 'intersects', aligned_array[j][0]

            substr = longest_substring(aligned_array[j][0],aligned_array[i][0])
            if len(substr) > 0:
                test2_failed = True
                print '*FAILED:', aligned_array[i][0], 'intersects', aligned_array[j][0]

    if len(aligned_array) > 1 and not test2_failed:
        if PRINT_PASS_AND_DATA:
            print 'PASSED: No common substrings'

def test_naive_multiple(rounds):
    rseq = ALPHABET

    for i in range(0,rounds):
        seq_reads = split_seq(rseq,10,5,6,1.0)
        test_naive_single(rseq,seq_reads)

test_naive_multiple(1000)


