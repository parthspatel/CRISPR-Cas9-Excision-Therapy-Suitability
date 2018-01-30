import pandas as pd

from Bio import SeqIO
from collections import Counter


class probe:

    dictDNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    degenerate = ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']

    def __init__(self, paths, length=23, frac=0.001):
        assert frac <= 1 and frac >= 0, 'Frac must be between 0 and 1 inclusively'
        try:
            for path in paths:
                f = open(path)
                f.close()
        except IOError:
            print('File is not accessible')

        self.paths = paths
        self.length = length
        self.frac = frac

    def importSequences(self):
        seqs = []
        for path in self.paths:
            with open(path) as handle:
                seqs += list(SeqIO.parse(handle, 'fasta'))
        self.seqs = seqs

        return self

    # Generate the kmer sequnces from the fasta input
    def generate(self):
        rmers = Counter()
        for seq in self.seqs:
            seqStr = str(seq.seq.ungap('-').upper())
            lenth = len(seqStr)

            for loc in range(len(seqStr) - self.length):
                rmers[seqStr[loc:loc + self.length]] += 1

        self.filterRare(rmers)
        self.filterPAM().filterDegenerate().convertPAM().extractProto()
        self.kmers = pd.DataFrame(self.kmers)

        return self

    # Remove Sequences that have low occurances (rare)
    def filterRare(self, rmers):
        cutoff = len(self.seqs) * self.frac
        nonRare = sum(count > cutoff for count in rmers.values())
        self.kmers = pd.Series([seq for seq, count in rmers.items() if count > cutoff])

        return self

    # Return the reverse compliment of the input sequence
    def reverseCompliment(self, seq):
        return ''.join([self.dictDNA[base] for base in list(seq[::-1])])

    # Convert any seqs starting in CC to their proper PAM seqs
    def convertPAM(self):
        self.kmers = [list(k) for k in self.kmers.values]
        for i, k in enumerate(self.kmers):
            l = ''.join(k)
            if l.startswith('CC'):
                self.kmers[i] = self.reverseCompliment(l)
            else:
                self.kmers[i] = l

        return self

    # Filter out any seqs that do no end or start with PAM or it's compliment
    def filterPAM(self):
        self.kmers = self.kmers[self.kmers.map(lambda x: x.endswith('GG') | x.startswith('CC'))]
        self.kmers.reset_index()
        return self

    # Filter out any seqs that contain degenerate bases
    def filterDegenerate(self):
        self.kmers = self.kmers[self.kmers.map(lambda x: all(b not in x for b in self.degenerate))]
        self.kmers.reset_index()
        return self

    # Split the protospacer from the PAM seqs
    def extractProto(self):
        pam = [0] * len(self.kmers)
        proto = [0] * len(self.kmers)
        for i, p in enumerate(self.kmers):
            l = ''.join(p)
            j = l[-3:]
            k = l[0: self.length - 3]

            pam[i] = j
            proto[i] = k

        self.pam = pd.DataFrame(pam)
        self.proto = pd.DataFrame(proto)

        return self

    # Print the input data to a csv file
    def toCSV(data, fileName):
        outFile = csv.writer(open(fileName, 'w'), delimiter=',', quoting=csv.QUOTE_ALL)
        for i in data:
            outFile.writerow(i)
        return
