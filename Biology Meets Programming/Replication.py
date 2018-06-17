#Find the reverse complement of a DNA string
def ReverseComplement(Pattern):
    revComp = ''
    for i in Pattern:
        if i == "A":
           revComp = "T" + revComp
        elif i == "T":
             revComp = "A" + revComp
        elif i == "G":
             revComp = "C" + revComp
        elif i == "C":
             revComp = "G" + revComp
    return revComp

def FrequentWords(Text, k):
    FrequentPatterns = [] # output variable
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    return FrequentPatterns

def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count
def PatternCount(Pattern, Text):
    count = 0 # output variable
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern :
            count = count + 1
            
    return count


def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern :
            positions.append(i)
    return positions


 def SymbolArray(Genome, symbol):
    array = {}
    # type your code here
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


 def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


 def Skew(genome):
    skew = {0:0}

    for i in range(1, len(genome)+1):

        if genome[i - 1] == "G":
            skew[i] = skew[i - 1] + 1
        elif genome[i - 1] == "C":
            skew[i] = skew[i - 1] - 1
        else:
            skew[i] = skew[i - 1]

    return skew


 def HammingDistance(p, q):
    hamming = 0
    for (i, x) in zip(p, q):
        if i != x:
           hamming += 1
    return hamming


 def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d :
           positions.append(i)
    return positions


 def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d :
           count += 1
    return count


# Finding most frequent k-mer with up to d mismatches
import collections

def AllKmers(Genome, k, d):
    kmers = []
    positions = []
    mismatches = []
    for i in range(len(Genome)- k + 1):
        kmers.append(Genome[i:i+k])
    for x in set (kmers):
        for y in kmers:
            if HammingDistance(x,y) <= d:
               positions.append(x)
    mismatch_count = collections.Counter(positions)
    for key,val in mismatch_count.items():
        if val == max(mismatch_count.values()):
           mismatches.append(key)
    return mismatches

