import random
import time
import math

random.seed()

letters_dna = ['A','T','G','C']

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
letters_alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
                    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def random_seq(letters,length):
    if length < 1:
        length = 1

    ret_str = ""

    for x in range(0,length):
        ret_str += letters[random.randint(0,len(letters)-1)]

    return ret_str

#Split a string into an array of random reads
#seq is the input string to split, length L
#num is the number of reads in the output, N
#smallest is the smallest read size, Rmin
#largest is the smallest read size, Rmax
#Accuracy will simulate errors in the read with rate A
def split_seq(seq, letters, num, smallest, largest, accuracy):
    if accuracy > 1.0 or accuracy < 0.5:
        return []

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
        new_read = seq[read_start_pos:read_end_pos+1]

        #introduce errors into reads
        if accuracy < 1.0:
            for i in range(0,len(new_read)):
                if random.randint(0,math.floor(1.0/(1.0-accuracy))) == 0:
                    nr_list = list(new_read)
                    nr_list[i] = letters[random.randint(0,len(letters)-1)]
                    new_read = "".join(nr_list)

        read_arr.append(new_read)

    return read_arr

#This algorithm runs in O(N^2 * L^2), where N=len(read_arr), L = max read size
def naive_polynomial_assemble(read_arr):
    begin_time = time.time()

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
        align_shift_arr = list()

        for r2 in range(r1+1,len(read_arr)):
            read2 = read_arr[r2]

            #Find the most probable alignment (largest # of matches) [O(L^2)]
            max_alignment = 0
            max_sum_of_squares = 0
            corresponding_shift = len(read1) #Nonsense value

            for shift in range(1-len(read2),len(read1)):
                alignment = 0
                sum_of_squares = 0
                sos_length = 0

                overlap1 = read1[max(0,shift):min(len(read1),shift+len(read2))]
                overlap2 = read2[max(0,-shift):min(len(read2),len(read1)-shift)]

                for i in range(0,len(overlap1)):
                    if overlap1[i] == overlap2[i]:
                        alignment += 1
                        sos_length += 1

                    if overlap1[i] != overlap2[i]:
                        sum_of_squares += sos_length*sos_length
                        sos_length = 0

                sum_of_squares += sos_length*sos_length

                if alignment > max_alignment:
                    max_alignment = alignment
                    #corresponding_shift = shift

                if sum_of_squares > max_sum_of_squares:
                    max_sum_of_squares = sum_of_squares
                    corresponding_shift = shift

                #Could add: if max_alignment = min(len(read1),len(read2)): break
                #Since it cannot be higher than the smaller of the lengths

            #If the max_alignment is insignificant, do not use the pair (r1,r2)
            #Otherwise, add it into the list of pairs to be merged
            align_shift_arr.append(corresponding_shift)

            #Determine the lower threshold for including into list
            #   This formula was essentially deduced empirically
            lower_threshold = math.pow(min(len(read1),len(read2)),1.5)
            if max_sum_of_squares > lower_threshold:
                align_score_arrs.append((max_sum_of_squares,r1,r2))

        align_shift_arrs.append(align_shift_arr)

    #Get alignments in the order they should be processed based on max score
    #This takes O(N^2), since the ao_first & ao_second have constant time 'in'
    #   and there are N^2 entries in align_score_arrs
    asa_sorted = sorted(align_score_arrs, key=lambda x: (x[0],-x[1],-x[2]))
    align_score_arrs = reversed(asa_sorted)
    align_score_arrs_copy = reversed(asa_sorted)

    align_order = list()
    ao_first = set()
    ao_second = set()

    for a in align_score_arrs:
        if a[0] == 0:
            if a[2] not in ao_first and a[2] not in ao_second:
                ao_second.add(a[2])
                align_order.append((len(read_arr),a[2]))
        else:
            ao_first.add(a[1])
            ao_second.add(a[2])

            align_order.append((a[1],a[2]))

    #Fill in shift array, which is triangular but needs to be square [0(N^2)]
    #Then align_shift_arrs[i][j] is the amount needed to shift j to align to i
    align_shift_arrs.append(list())

    for j in range(0,len(align_shift_arrs)):
        align_shift_arrs[j].insert(0,0)
        for i in range(0,j)[::-1]:
            align_shift_arrs[j].insert(0,-align_shift_arrs[i][j-len(read_arr)])

    if VERBOSE:
        print 'Alignment-score array:'
        for a in align_score_arrs_copy:
            print a
        print ''

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
    #I believe this runs in O(NlogN)
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
                                insert_into_mergelet(mergelet_arr[i],
                                                     insert_read[0],
                                                     insert_read[1]+shift1-shift2)

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
    #Runs in O(N*L)
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

    print "Total time: " + str(time.time() - begin_time) + "s"

    return aligned_array

#Get the sum-squared score for a string, roughly representing its entropy
#1.0 means no sequential sequences, >>1.0 means long sequential sequences
def get_sum_squared_score(instr):
    char_last = ''
    count = 0
    sos_count = 0
    for i in range(0,len(instr)):
        if instr[i] == char_last:
            count += 1
        else:
            sos_count += count*count
            count = 1
            char_last = instr[i]
    
    sos_count += count*count
    return float(sos_count)/(len(instr))

def test_naive_single(rseq,seq_reads):
    PRINT_DATA = True

    if PRINT_DATA:
        print 'seq:',rseq
        #print 'data:',seq_reads,'\n'

        avg_length = 0
        for lr in seq_reads:
            avg_length += len(lr)
        avg_length /= len(seq_reads)

        print 'coverage: ' + str(math.floor(avg_length*len(seq_reads)/len(rseq))) + 'x'

    #Get assembled string
    aligned_array = naive_polynomial_assemble(seq_reads)

    #Check assembled contigs against original sequence
    rseq_coverage = 0
    for i in range(0,len(aligned_array)):
        read1 = aligned_array[i][0]
        read2 = rseq

        #Find the most probable alignment (largest # of matches) between read1&2
        max_alignment = 0
        corresponding_shift = len(read1) #Nonsense value

        for shift in range(1-len(read2),len(read1)):
            alignment = 0

            overlap1 = read1[max(0,shift):min(len(read1),shift+len(read2))]
            overlap2 = read2[max(0,-shift):min(len(read2),len(read1)-shift)]

            for j in range(0,len(overlap1)):
                if overlap1[j] == overlap2[j]:
                    alignment += 1

            if alignment > max_alignment:
                max_alignment = alignment

        rseq_coverage += max_alignment

        if PRINT_DATA:
            perc_match = math.floor(100*max_alignment/len(read1))
            print "read",str(i)+"(acc="+str(perc_match)+"%):",read1

    if PRINT_DATA:
        rseq_match = math.floor(100*rseq_coverage/len(rseq))
        print "Total  acc="+str(rseq_match)+"%"
        print "SOS score:",get_sum_squared_score(rseq)
    return aligned_array

def test_naive_multiple(rounds):

    for i in range(0,rounds):
        seq = random_seq(letters_dna,100)
        out = test_naive_single(seq,split_seq(seq,letters_dna,300,10,10,1.0))

VERBOSE = False #For printing intermediate steps and debugging
test_naive_multiple(1)

