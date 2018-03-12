peakheaps
=========

![peakheapsABC.png](peakheapsABC.png)

Example:

    $ ./peakheaps.py condA.bed condB.bed condC.bed > peakheapsABC.bed
    2 peaks found in condA.bed
    2 peaks found in condB.bed
    2 peaks found in condC.bed
    3 connected components 
    
    $ cat peakheapsABC.bed 
    chrI	100	200	site1_condA,site1_condB,site1_condC	3	.
    chrI	255	355	site3_condA,site3_condB	2	.
    chrI	240	270	site2_condC	1	.

(So output coordinates are calculated as the mean of the constituent peaks, name contains comma-separated names of the constituent peaks, score (=5th column) is the size of the connected component (=number of constituent peaks).