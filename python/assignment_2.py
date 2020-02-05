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



from python.assignment_1 import read_genome, naive_2mm



'''
Implement the pigeonhole principle using Index class
to find exact matches for the partitions. Assume P always has 
length 24, and that we are looking for approximate matches with
up to 2 mismatches (substitutions). We will use an 8-mer index.

Download the Python module for building a k-mer index.

https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py

Write a function that, given a length-24 pattern P and given an 
Index\verb|Index|Index object built on 8-mers, finds all approximate 
occurrences of P within T with up to 2 mismatches. Insertions and 
deletions are not allowed. Don't consider any reverse complements.

How many times does the string GGCGCGGTGGCTCACGCCTGTAAT\verb|
GGCGCGGTGGCTCACGCCTGTAAT|GGCGCGGTGGCTCACGCCTGTAAT, which is 
derived from a human Alu sequence, occur with up to 2 substitutions 
in the excerpt of human chromosome 1? (Don't consider reverse complements here.)

Hint 1: Multiple index hits might direct you to the same match multiple 
times, but be careful not to count a match more than once.

Hint 2: You can check your work by comparing the output 
of your new function to that of the naive_2mm function 
implemented in the previous module.
'''

import bisect

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits



def build_index_and_call_pattern(p, t):
  
  # create an index of 8-mers for the text
  
  index = Index(t, 8)
  
  # now for each 8 mer in a length 24 pattern we will query the data
  # we are allowing up to 2 mismatches, so we must match at least one
  # of the 3 8-mer partitions in the 24 bp pattern
  
  # for testing
  from python.assignment_1 import naive_2mm
  hit_idxs = []
  verified_hit_idxs = []
  
  for i in range(0, 24-8+1, 8): 
    kmer = p[i:i+8]
    print('kmer: ', kmer)
    res = index.query(kmer)
    print(res)
    if len(res) >= 1: 
      hit_idxs.extend(res)
      for match_ind in res: 
        # verify match
        ind = match_ind-i
        if len(naive_2mm(p, t[ind:ind+24])) >= 1:
          verified_hit_idxs.append(ind)
  
  print('Total kmer hits: ', nhits)
  print('Total amount of unique indices returned: ', len(set(hit_idxs)))
  print('Unique verified indexed hit counts: ', len(set(verified_hit_idxs)))
  
  print('Naive hit counts: ', len(naive_2mm(p, t)))
  
t = read_genome('data/week_2_hw.fasta')
build_index_and_call_pattern('GGCGCGGTGGCTCACGCCTGTAAT', t)


import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits




'''
Write a function that, given a length-24 pattern P and 
given a SubseqIndex object built with k = 8 and ival = 3, 
finds all approximate occurrences of P within T 
with up to 2 mismatches.

When using this function, how many total index hits are there when 
searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions 
in the excerpt of human chromosome 1? (Again, don't consider 
reverse complements.
'''


def build_index_and_call_pattern_2(p, t): 

  index = SubseqIndex(t, 8, 3)
  
  hit_idxs = []
  verified_hit_idxs = []
  
  for i in range(0, 3):
    kmer = p[i:]
    res = index.query(kmer)
    print(res)
    if len(res) >= 1:
      hit_idxs.extend(res)
      for match_ind in res:
        # verify match
        ind = match_ind-i
        print(t[ind:ind+24])
        if len(naive_2mm(p, t[ind:ind+24])) >= 1:
          verified_hit_idxs.append(ind)
  
  
  
  print('Number of index hits: ', len(set(hit_idxs)))
  print('Verified indexed hit counts: ', len(set(verified_hit_idxs)))
  print('Naive hit counts: ', len(naive_2mm(p, t)))


p = 'GGCGCGGTGGCTCACGCCTGTAAT'
t = read_genome('data/week_2_hw.fasta')

build_index_and_call_pattern_2(p, t)

'''
messed up number of hits, was counting number of indices
I was checking but there could be multiple indices per hit
'''




