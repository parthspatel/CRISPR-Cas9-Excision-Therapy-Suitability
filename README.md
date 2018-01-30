# Senior Design: CRISPR-Cas9 Excision Therapy Suitability


This is designed to determine if a patient is suitable for a *Sp*Cas9 therapy.

First, it finds all the 23bp sequences within the Los Alamos Database. Once the algorithm has picked out the necessary sequences the estimated number of probes comes out to 44,772. However, this number can be further reduced by filtering out rare sequences, sequences that only occur within the patient database less than 0.1% of the time. These sequences are not needed in the final micro array because the data they provide represents an extremely small portion of the total set and are not essential to final analysis of the array. By filtering out the rare sequences the number of probes goes from 44,772 to 4,115. Being able to reduce the number of probes on the micro array allows for a reduced cost. And since having fewer probes reduces the amount of space each sample takes up on the array, multiple samples can be tested at once. This may allow for replicate testing of the same patient sample or multiple patient testing done on one array.

Next, it will use the crisprtree module to emulate *Sp*Cas9 binding and determine which probes can bind to which DNA sequences.  In order to determine patient suitability this information will be combinded with a microarray analysis.  This analysis is still currently in development.

To learn more about the project, please refrence my webpage: [parthspatel.me](http://parthspatel.me/senior-design/)
The crisprtree python module has been developed by Dr. Will Dampier, to refrence my fork: [crisprtree](https://github.com/parthspatel/crisprtree)
