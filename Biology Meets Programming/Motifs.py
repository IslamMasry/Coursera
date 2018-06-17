def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for symbol in count:
        profile[symbol] = count[symbol]
    for i in "ACGT":
        for j in range(k):
            profile[i][j] /= t
    return profile  



def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    score = 0
    t = 0
    for j in range(k):
        m = 0
        for symbol in "ACGT":                                  
            score += count[symbol][j]
            if count[symbol][j] > m: 
               m = count[symbol][j]
        score -= m
    return score



def Pr(Text, Profile):
    p = 1
    for i in range (len(Text)):
        p *= Profile[Text[i]][i]
    return p


def ProfileMostProbablePattern(Text, k, Profile):
    pattern = []
    m = 0
    maxx = ""
    for i in range(len(Text)-k+1):
        pattern.append(Text[i:i+k])
        if Pr(pattern[i], Profile) > m:
           m = Pr(pattern[i], Profile)
           maxx = pattern[i]
    if m <= 0:
       maxx = pattern[0]
    return maxx


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
          BestMotifs = Motifs
    return BestMotifs


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol in count:
        profile[symbol] = count[symbol]
    for i in "ACGT":
        for j in range (k):
            profile[i][j]/= t+4
        
    return profile


def Motifs(Profile, Dna):
    most = []
    for i in range(0, len(Dna)):
        most.append(ProfileMostProbablePattern(Dna[i], Profile))
    return most


def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
          return BestMotifs


def RandomMotifs(Dna, k, t):
    r = []
    for i in range (0, len(Dna)):
        s = Dna[i]
        x = randint(0, len(Dna[0])-k)
        r.append(s[x:x+k])
    return r

def RepeatedRandomizedMotifSearch(Dna, k, t):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(1000):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs


def WeightedDie(Probabilities):
    kmer = '' # output variable
    x = {}
    t = 0.0
    P = random.uniform(0, 1)
    for key, value in Probabilities.items():
        x[key] = []

    for key, value in Probabilities.items():
        x[key].append(t)
        x[key].append(t+Probabilities[key])
        t += Probabilities[key]
    for key, value in x.items():
        if P > x[key][0] and P <= x[key][1]:
           kmer = key
        
    return kmer


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for j in range(1, N):
        i = randint(0, t)
        profile = ProfileWithPseudocounts(M)
        
        if Score(M) < Score(BestMotifs):
           BestMotifs = M
    return BestMotifs

