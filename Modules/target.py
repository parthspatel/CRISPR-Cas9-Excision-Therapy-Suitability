import pandas as pd

from probe import probe
from Bio import SeqIO
from collections import Counter


class target:

    dictDNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}
    degenerate = ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']

    def __init__(self, data, paths, length=23):
        try:
            for path in paths:
                f = open(path)
                f.close()
        except IOError:
            print('File is not accessible')

        self.paths = paths
        self.length = length

        if type(data) is probe:
            self.size = len(data.proto)
        elif type(data) is pd.core.frame.DataFrame:
            self.size = len(data)
        else:
            self.size = data

    # Generate the target sequences
    def generate(self):
        seqs = []
        LTRBeg = None
        LTREnd = None
        for path in self.paths:
            with open(path) as handle:
                seqs += list(SeqIO.parse(handle, 'fasta'))
            if not LTRBeg:
                LTRBeg = seqs[-1]
            elif not LTREnd:
                LTREnd = seqs[-1]
        self.LTR = [LTRBeg, LTREnd]

        begKmers = self.convertPAM(self.filterDegenerate(
            self.filterPAM(self.generateKmers(self.LTR[0]))))
        endKmers = self.convertPAM(self.filterDegenerate(
            self.filterPAM(self.generateKmers(self.LTR[1]))))

        self.kmers = [begKmers, endKmers]
        return self

    # Generate the kmers from the LTR sequnces
    def generateKmers(self, seq):
        rmers = Counter()
        seqStr = str(seq.seq)  # .ungap('-').upper()
        for loc in range(len(seqStr) - self.length):
            rmers[seqStr[loc:loc + self.length]] += 1
        # kmers = []
        # for mer in rmers:
        #     if '-' not in mer:
        #         kmers += [str(mer)]
        return pd.Series(rmers)

    # Filter out any seqs that do no end or start with PAM or it's compliment
    def filterPAM(self, kmers):
        kmers = pd.Series(kmers.index)
        kmers = kmers[kmers.map(lambda x: x.endswith('GG') | x.startswith('CC'))]
        kmers.reset_index()
        return kmers

    # Filter out any seqs that contain degenerate bases
    def filterDegenerate(self, kmers):
        kmers = kmers[kmers.map(lambda x: all(b not in x for b in self.degenerate))]
        kmers.reset_index()
        return kmers

    # Convert any seqs starting in CC to their proper PAM seqs
    def convertPAM(self, kmers):
        kmers = [list(k) for k in kmers.values]
        for i, k in enumerate(kmers):
            l = ''.join(k)
            if l.startswith('CC'):
                kmers[i] = self.reverseCompliment(l)
            else:
                kmers[i] = l
        return kmers

    # Return the reverse compliment of the input sequence
    def reverseCompliment(self, seq):
        return ''.join([self.dictDNA[base] for base in list(seq[::-1])])

    # Print the input data to a csv file
    def toCSV(data, fileName):
        outFile = csv.writer(open(fileName, 'w'), delimiter=',', quoting=csv.QUOTE_ALL)
        for i in data:
            outFile.writerow(i)
        return
