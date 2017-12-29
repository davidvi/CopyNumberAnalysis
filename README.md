# CopyNumberAnalysis
seqCNA C++ code to generate files for readSeqsumm

The R project can be found here:
[Bioconductor: seqCNA](https://bioconductor.org/packages/release/bioc/html/seqCNA.html)

This is the C++ code from the R package David Mosen-Ansorena published, adapted to run from to command line. I used this to generate the files used for down-steam analysis of copy number variations using the seqCNA R package as the runSeqsumm function did not work on my machine.

### to build and run the binary
```
git clone https://github.com/davidvi/CopyNumberAnalysis.git
cd CopyNumberAnalysis
g++ copy_number_analysis.cpp -o cna
./cna -w 50 -i input.bam -o output.txt
```

options:
-w window size (*1000)
-i input bam file
-o output txt file for use in seqCNA R package
