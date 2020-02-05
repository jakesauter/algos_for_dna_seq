def read_genome(filename):
  genome = ''
  with open(filename, 'r') as f:
    for line in f:
      # ignore header line with genome information
      if not line[0] == '>':
        genome += line.rstrip()
  
  return genome
    
def read_fastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive(p, t):
  occurrences = []
  for i in range(len(t) - len(p) + 1):
    # loop over alignments
    match = True
    for j in range(len(p)):
      # loop over characters
      if t[i + j] != p[j]:
        # compare characters
        match = False
        break
    if match:
        occurrences.append(i)  # all chars matched; record
  
  return occurrences
    

def naive_with_rc(p, t):
  occurrences = naive(p,t)
  occurrences.extend(naive(reverse_complement(p), t))
  return list(sorted(set(occurrences)))


def naive_2mm(p, t): 
  occurrences = []
  # for all valid alignments
  for i in range(len(t) - len(p) + 1):
    match = True
    mm = 0
    # for all characters in pattern
    for j in range(len(p)):
      if t[i + j] != p[j]:
        # compare characters
        mm += 1
        if mm > 2: 
          match = False
          break
    if match:
        occurrences.append(i)  # all chars matched; record
  
  return occurrences
    

# data = read_genome('data/lambda_virus.fa')
# 
# len(naive_with_rc('AGGT', data))
# len(naive_with_rc('TTAA', data))
#   
# min(naive_with_rc('ACTAAGT', data))
# min(naive_with_rc('AGTCGA', data))
#   


# data, quality = read_fastq('data/first1000.fastq')
#   
# bad_reads_per_cycle = []
# 
# for i in range(100):
#   bad_reads_per_cycle.append(len([q[i] for q in quality if q[i] == '#']))
# 
#   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

