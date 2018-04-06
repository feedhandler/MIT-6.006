#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *
from collections import defaultdict
from collections import deque

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
      self.d = defaultdict(list)
      for k,v in pairs:
        self.put(k,v)
        
    # Associates the value v with the key k.
    def put(self, k, v):
      self.d[k].append(v)
        
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
      return self.d[k]
      
    # Try to print it
    def __str__(self):
      return str(self.d)
        
# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
  subseq = deque()
  rh = None # this will be our rolling hash
  pos = -1 # zero-indexed position of the subsequence

  for newItem in seq:
    subseq.append(newItem) # append to right side
    if len(subseq)<k:
      continue
      
    pos += 1
      
    if rh is None:
      rh = RollingHash(subseq)
    else:
      oldItem = subseq.popleft()
      rh.slide(oldItem, newItem)

    yield (rh.current_hash(), deque(subseq), pos)
    
# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
  i=0
  # TODO this is a bit inefficient because we are computing the
  # hash for every subsequence but then we're only yielding
  # every m-th one. So we are wasting quite a lot of effort
  for hsh,subseq,pos in subsequenceHashes(seq, k):
    if i % m == 0:
      yield hsh,subseq,pos
    i += 1

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
# Return pairs of offsets into the inputs
# A tuple(x, y) being re-turned indicates that the k-length subsequence at position
# x in the first input matches the subsequence at position y in the second input
def getExactSubmatches(a, b, k, m):

  mdA = Multidict()
  for hsh,subseq,pos in intervalSubsequenceHashes(a, k, m):
    mdA.put(hsh, (subseq, pos))
    
  for hashB,subseqB,posB in subsequenceHashes(b, k):
    for subseqA, posA in mdA.get(hashB):
      if subseqA == subseqB:
        yield posA, posB
        

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0])
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
