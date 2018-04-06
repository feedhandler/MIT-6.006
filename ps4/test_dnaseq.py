from dnaseq import *

### Testing ###

class TestSubsequenceHashes(unittest.TestCase):
  def test_subsequenceHashes(self):
    seq=['A','9','8','7','6','5','4','3','2','1']
    expected = [
      (178799, deque('A9876'), 0),
      (159191, deque('98765'), 1),
      (156390, deque('87654'), 2),
      (153589, deque('76543'), 3),
      (150788, deque('65432'), 4),
      (147987, deque('54321'), 5),
      ]
    expGen = (e for e in expected)
    for x in subsequenceHashes(iter(seq),5):
      self.assertEqual(next(expGen), x)
      
#~ class TestIntervalSubsequenceHashes(unittest.TestCase):
  #~ def test_intervalSubsequenceHashes(self):
    #~ seq=['F','E','D','C','B','A','9','8','7','6','5','4','3','2','1','0']
    #~ for x in intervalSubsequenceHashes(iter(seq),5,3):
      #~ print x
              
class TestRollingHash(unittest.TestCase):
    def test_rolling(self):
        rh1 = RollingHash('CTAGC')
        rh2 = RollingHash('TAGCG')
        rh3 = RollingHash('AGCGT')
        rh1.slide('C','G')
        self.assertTrue(rh1.current_hash() == rh2.current_hash())
        rh1.slide('T','T')
        self.assertTrue(rh1.current_hash() == rh3.current_hash())

class TestMultidict(unittest.TestCase):
    def test_multi(self):
        foo = Multidict()
        foo.put(1, 'a')
        foo.put(2, 'b')
        foo.put(1, 'c')
        self.assertTrue(foo.get(1) == ['a','c'])
        self.assertTrue(foo.get(2) == ['b'])
        self.assertTrue(foo.get(3) == [])

# This test case may break once you add the argument m (skipping).
class TestExactSubmatches(unittest.TestCase):
   def test_one(self):
       foo = 'yabcabcabcz'
       bar = 'xxabcxxxx'
       matches = list(getExactSubmatches(iter(foo), iter(bar), 3, 1))
       correct = [(1,2), (4,2), (7,2)]
       self.assertTrue(len(matches) == len(correct))
       for x in correct:
           self.assertTrue(x in matches)

unittest.main()
