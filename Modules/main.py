import numpy as np
import pandas as pd
import csv

from probe import probe
from target import target
from crisprtree import preprocessing
from crisprtree import estimators
from crisprtree import evaluators

# File Paths to the input data
# Use argparse to dynamically input paths
paths = ['C:\\Users\\parth\\Desktop\\Python Workspace\\Senior Design\\Data\\hiv-1-700.fixed.fst',
         'C:\\Users\\parth\\Desktop\\Python Workspace\\Senior Design\\Data\\hiv-9086-9717.fixed.fst']

# Generate the probes/kmers for the Microarray & seqs that represent potential seqs for CRISPR
p = probe(paths).importSequences().generate()

# Print the probes to a csv file
# p.toCSV(p.kmers[0],'C:\Users\parth\Desktop\Python Workspace\Senior Design\Data\probes.csv')


# Generate the target sequences from the refrences HXB2 cell line
t = target(p.proto, paths).generate()


# Filter out target seqs that have missing bases '-' and create all possible valid pairs between the protospacer and the target sequnces.
#-------------------------------------------------------------------------------
inputSeqs = pd.DataFrame()
#p.proto = p.proto[0:1000]
for i, spacer in enumerate(p.proto[0]):
    if i % 1000 == 0:
        print(spacer, i)
    # Runs on the begining of the LTR (t.kmers[0]), the end of the LTR (t.kmers[1])
    for _, mer in enumerate(t.kmers[1]):
        if '-' not in mer:
            inputSeqs = inputSeqs.append(pd.DataFrame([spacer], [mer]))

inputSeqs = inputSeqs.reset_index()
inputSeqs.columns = ['target', 'gRNA']
inputSeqs = inputSeqs.reindex(columns=['gRNA', 'target'])


# Run the pairs of sequence through the binding estimator
#-------------------------------------------------------------------------------

# Put the input data in the format required by crisprtree
pipe = np.array(inputSeqs)

# Test using the MIT estimator, upgrade to CFD in the future
est = estimators.MITEstimator(pipe)
# to get the numerical scores use predict_proba, for boolean use predict
estResults = est.build_pipeline().predict(pipe)

pipe = pd.DataFrame(pipe)

# Append the results from the estimator and reset the data frame column titles
pipe['binding'] = list(estResults)
pipe.columns = ['gRNA', 'target', 'binding']

# Filter out and print any sequences that can bind
bindable = pipe[pipe['binding'] == True].reset_index()

print(bindable)


# Results:
# Total bindable probes: 126 probes + all HXB2
#-------------------------------------------------------------------------------


# Begining of the LTR:
#     index                  gRNA                   target  binding
# 0     590  CAGCTGCTTTTTGCCTGTAC  CAGCTGCTTTTTGCCTGTACTGG     True
# 1     602  AGCTGCTTTTTGCCTGTACT  AGCTGCTTTTTGCCTGTACTGGG     True
# 2     629  TTGCCTGTACTGGGTCTCTC  TTGCCTGTACTGGGTCTCTCTGG     True
# 3     637  TAACCAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 4     655  GGTTAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True
# 5     670  GTTAGACCAGATCTGAGCCT  GTTAGACCAGATCTGAGCCTGGG     True
# 6     677  AGAGCTCCCAGGCTCAGATC  AGAGCTCCCAGGCTCAGATCTGG     True
# 7     687  ATCTGAGCCTGGGAGCTCTC  ATCTGAGCCTGGGAGCTCTCTGG     True
# 8    1247  ATTTGAGCCTGGGAGCTCTC  ATCTGAGCCTGGGAGCTCTCTGG     True
# 9    1335  AGAACTCCCAGGCTCAGATC  AGAGCTCCCAGGCTCAGATCTGG     True
# 10   1662  GACAAGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 11   2632  AGCTGCTTTCTGCCTGTACT  AGCTGCTTTTTGCCTGTACTGGG     True
# 12   2653  TAACAAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 13   2671  TGTTAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True
# 14   3283  TAGCTAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 15   3595  AGCTAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True
# 16   3610  GCTAGACCAGATCTGAGCCT  GTTAGACCAGATCTGAGCCTGGG     True
# 17   3643  CAAGGATATCTTGTCTTCGT  CAAGGATATCTTGTCTTCGTTGG     True
# 18   3658  CAGGGAAGTAGCCTTGTGTG  CAGGGAAGTAGCCTTGTGTGTGG     True
# 19   4004  AGCTGCTCTTTGCCTGTACT  AGCTGCTTTTTGCCTGTACTGGG     True
# 20   4790  CAGCCGCTTTTTGCCTGTAC  CAGCTGCTTTTTGCCTGTACTGG     True
# 21   4802  AGCCGCTTTTTGCCTGTACT  AGCTGCTTTTTGCCTGTACTGGG     True
# 22   5424  AAGGATATCTTGTCTTCGTT  AAGGATATCTTGTCTTCGTTGGG     True
# 23   5440  TGTGGTAGATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 24   5455  CTGTGGATCTACCACACACA  CTGTGGATCTACCACACACAAGG     True
# 25   6660  GACAAGACATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 26   6672  TGTTGTAGATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 27   7153  TCGCTTGTACTGGGTCTCTC  TTGCCTGTACTGGGTCTCTCTGG     True
# 28   7265  TCGCCTGTACTGGGTCTCTC  TTGCCTGTACTGGGTCTCTCTGG     True
# 29   7528  GAAAAGAGATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# ..    ...                   ...                      ...      ...
# 69  28764  TGTGATAGATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 70  29424  GGCAAGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 71  29586  CTGGGAAGTAGCCTTGTGTG  CAGGGAAGTAGCCTTGTGTGTGG     True
# 72  30006  CAGGGAAATAGCCTTGTGTG  CAGGGAAGTAGCCTTGTGTGTGG     True
# 73  30979  AGGTAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True
# 74  31424  TATGGTAGACCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 75  32085  GGCTAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True
# 76  32162  CGGGGAAGTAGCCTTGTGTG  CAGGGAAGTAGCCTTGTGTGTGG     True
# 77  37599  CTATGGATCTACCACACACA  CTGTGGATCTACCACACACAAGG     True
# 78  38059  TAGCCAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 79  38561  AGAGCTCTCAGGCTCAGATC  AGAGCTCCCAGGCTCAGATCTGG     True
# 80  39250  TGTGGTAAATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 81  39823  TACCTAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 82  40538  TGTTATAGACCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 83  40624  GACAGGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 84  41478  GACGAGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 85  42738  GACGAGAGATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 86  43044  TGTGGTACATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 87  43228  GACAGGAGATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 88  43284  GAAGAGAGATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 89  44110  GAAGAGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 90  44739  ATGTGGATCTACCACACACA  CTGTGGATCTACCACACACAAGG     True
# 91  45788  TGTTGTACATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 92  45970  TGTGGTACACCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 93  47875  TTGTGGATCCACCACACACA  CTGTGGATCTACCACACACAAGG     True
# 94  48125  TGACTAGAGAGACCCAGTAC  TAACCAGAGAGACCCAGTACAGG     True
# 95  48294  TGTTATAGATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 96  48350  TGTGATACATCCACAGATCA  TGTGGTAGATCCACAGATCAAGG     True
# 97  50214  GACGGGATATCCTTGATCTG  GACAAGATATCCTTGATCTGTGG     True
# 98  51069  AGTCAGACCAGATCTGAGCC  GGTTAGACCAGATCTGAGCCTGG     True


#-------------------------------------------------------------------------------


# End of the LTR:
#     index                  gRNA                   target  binding
# 0      25  AAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 1      33  TAGTTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 2     395  TAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 3     403  TAGCTTGTAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 4     655  AAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 5     673  TAGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 6    1436  CACTGACTAAAAGGGTCTGA  CACTGACTAAAAGGGTCTGAGGG     True
# 7    1442  ACACTGACTAAAAGGGTCTG  ACACTGACTAAAAGGGTCTGAGG     True
# 8    1644  TCAGACCCTTTTAGTCAGTG  TCAGACCCTTTTAGTCAGTGTGG     True
# 9    2569  TCAGACCATTTTAGTCAGTG  TCAGACCCTTTTAGTCAGTGTGG     True
# 10   3333  TAACTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 11   4523  CAACTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 12   5773  CAGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 13   7368  TAGTTTGTAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 14   7568  TAACTTGTAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 15   8413  AAGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 16   8435  AAACACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 17   8443  TAGCTTGAAACACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 18   8688  TAGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 19   9700  TAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 20  11525  AAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 21  11821  CGCTAACTAAAAGGGTCTGA  CACTGACTAAAAGGGTCTGAGGG     True
# 22  11949  TCAGACCTTTTTAGTCAGTG  TCAGACCCTTTTAGTCAGTGTGG     True
# 23  12171  CACTAACTAAAAGGGTCTGA  CACTGACTAAAAGGGTCTGAGGG     True
# 24  12249  TCAGACCATCTTAGTCAGTG  TCAGACCCTTTTAGTCAGTGTGG     True
# 25  17173  TTGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 26  17238  GAGCTTGAAGCACCATCCAA  TAGCTTGTAGCACCATCCAAAGG     True
# 27  18540  CAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
# 28  19000  CAGCACCATCCAAAGGTCAG  TAGCACCATCCAAAGGTCAGTGG     True
