'''
Implement versions of the naive exact matching and Boyer-Moore 
algorithms that additionally count and return (a) the number of character
comparisons performed and (b) the number of alignments tried. Roughly speaking, 
these measure how much work the two different algorithms are doing.

For a few examples to help you test if your enhanced versions of the naive exact 
matching and Boyer-Moore algorithms are working properly, see these notebooks:


from naive_with_counts import naive_with_counts
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

([40], 41, 46)

p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

([0, 19], 20, 35)


# Implement boyer_moore_with_counts by extending boyer_moore function
from bm_with_counts import boyer_moore_with_counts
from bm_preproc import BoyerMoore

p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)


# character comparisons high: [40] 12 11
([40], 12, 15)

p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)

# character comparison low: [0, 19] 5 3
([0, 19], 5, 18)

'''


'''
Notes: 

So it seems that we do not need to implement anything new here, just need to modify the algorithms
to determine 
 a) number of character comparisons performed 
 b) the number of alignments tried
'''


from python.boyer_moore_pre_proc import BoyerMoore

def naive_with_counts(p, t):
  character_comparisons = 0
  alignments_tried = 0
  
  occurrences = []
  for i in range(len(t) - len(p) + 1):
    # loop over alignments
    alignments_tried += 1
    match = True
    for j in range(len(p)):
      # loop over characters
      character_comparisons+=1
      if t[i + j] != p[j]:
        # compare characters
        match = False
        break
    if match:
        occurrences.append(i)  # all chars matched; record
  
  return (occurrences, alignments_tried, character_comparisons)

def boyer_moore_with_counts(p, p_bm, t):
    '''
    Do Boyer-Moore matching. p=pattern, t=text,
    p_bm=BoyerMoore object for p 
    '''
    character_comparisons = 0
    alignments_tried = 0
    
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        alignments_tried += 1
        for j in range(len(p)-1, -1, -1):
          character_comparisons += 1
          if p[j] != t[i+j]:
            skip_bc = p_bm.bad_character_rule(j, t[i+j])
            skip_gs = p_bm.good_suffix_rule(j)
            shift = max(shift, skip_bc, skip_gs)
            mismatched = True
            break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return (occurrences, alignments_tried, character_comparisons)



from python.assignment_1 import read_fastq


t = read_fastq('data/week_2_hw.fasta')
