import random

random.seed()

VERBOSE = True #For printing intermediate steps
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

    #Merge mergelets if possible, maintaining sorted order through a linear walk

    #Create list of mergelets to merge, and a separate list of pairs and shifts
    merge_mergelets = list()
    for i in range(0,len(mergelet_arr)):
        for read1 in mergelet_arr[i]:
            for j in range(i+1,len(mergelet_arr)):
                for read2 in mergelet_arr[j]:
                    if read1[0] == read2[0]:
                        if not [i,j,read1[1]-read2[1]] in merge_mergelets:
                            merge_mergelets.append([i,j,read1[1]-read2[1]])

    #To maintain topological order
    merge_mergelets = merge_mergelets[::-1]

    if VERBOSE:
        print("Merging mergelets:")

    #These lists will be used for accurately collapsing merges
    merge_mergelet_pairs = list()
    merge_mergelet_shifts = list()

    #Go through mergelets one at a time, merging them and collapsing the others
    #   to maintain consistency.
    while len(merge_mergelets) > 0:
        if VERBOSE:
            print 'To merge: ' + str(merge_mergelets)

        #Get first mergelet in list
        (ml1,ml2,all_shift) = merge_mergelets.pop(0)

        #Log this, in case other merges use ml1 or ml2
        merge_mergelet_pairs.append((ml1,ml2))
        merge_mergelet_shifts.append(all_shift)
        merge_mergelet_pairs.append((ml2,ml1))
        merge_mergelet_shifts.append(-all_shift)

        for (read,shift) in mergelet_arr[ml2]:
            if not (read, shift+all_shift) in mergelet_arr[ml1]:
                insert_into_mergelet(mergelet_arr[ml1], read, shift+all_shift)

        mergelet_arr.remove(mergelet_arr[ml2])

        if VERBOSE:
            print 'Merged ' + str([ml1,ml2,all_shift]) + ': ' + str(mergelet_arr)

        #Collapse merges:
        #After the current merge, all merge requests will be subsequently
        #   shifted down by 1, and the shifts changed to reflect the prior merge
        for i in range(0,len(merge_mergelets)):
            if merge_mergelets[i][0] >= ml2:
                pair_change = (merge_mergelets[i][0], merge_mergelets[i][0]-1)
                if pair_change in merge_mergelet_pairs:
                    consistency_shift = merge_mergelet_shifts[merge_mergelet_pairs.index(pair_change)]
                    merge_mergelets[i][2] += consistency_shift

                merge_mergelets[i][0] -= 1

            if merge_mergelets[i][1] >= ml2:
                pair_change = (merge_mergelets[i][1], merge_mergelets[i][1]-1)
                if pair_change in merge_mergelet_pairs:
                    consistency_shift = merge_mergelet_shifts[merge_mergelet_pairs.index(pair_change)]
                    merge_mergelets[i][2] += consistency_shift

                merge_mergelets[i][1] -= 1

            #Mark a redundant merge request for deletion
            if merge_mergelets[i][0] == merge_mergelets[i][1]:
                merge_mergelets[i] = [0,0,0]

        #Delete finished merge requests
        while [0,0,0] in merge_mergelets:
            merge_mergelets.remove([0,0,0])

    if VERBOSE:
        print ''
        print 'Final mergelet array: '
        print str(mergelet_arr)
        print ''

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
        aligned_array.append(list())
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

def main():
    rseq = ALPHABET
    seq_reads = split_seq(rseq,10,5,10,1.0)
    naive_polynomial_align(seq_reads)

main()
