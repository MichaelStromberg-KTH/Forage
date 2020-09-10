![Forage](https://github.com/MichaelStromberg-KTH/Forage/raw/main/forage.png)

Back in 2004, I created a Bayesian variant caller that used two different neural networks to filter false positives. The first NN uses the LVQ3 algorithm and the second NN uses the Optimal Brain Damage (OBD) algorithm. That algorithm name brings a smile to my face every time I say it.

These were very different times - long before BAM files and VCF files existed. Instead, we often used phred and ace files for the reads and alignments. This tool produced the output directly to stdout.

For more esoteric trivia, previous versions of Forage worked directly on assemblies produced by the Paracel EST assembler.

## Demo

Just run `make` and `./RunDemo.sh`

The resulting output looks like this:

```
Forage: SNP Discovery Software (c) 2004 Michael Stromberg
=========================================================
Looking for paralogs in Contig1
  - zb80g02.s1 marked as a paralog: P(NAT): 0.079789
  - zc43a10.s1 marked as a paralog: P(NAT): 0.184623
  - zn39a09.s1 marked as a paralog: P(NAT): 0.749026
Looking for SNPs in Contig1 (28 sequences)
Forage SNP @ 38 (unpadded: 38) P(SNP): 1.0000 P(VAR): 1.0000
       Freq3: 0.79 | B3: 49 B2: 37 | A: 4 C: 0 G: 0 T: 15
       CAACCATTTTTTTTTTTTTTTTTTTTTTTTTT W TAAAACAGTAGAAACAAGGTTGACTTTATTCC

Forage SNP @ 512 (unpadded: 469) P(SNP): 0.9909 P(VAR): 0.9909
       Freq3: 0.80 | B3: 34 B2: 45 | A: 3 C: 0 G: 12 T: 0
       CTAGAGGTTAAAAATGACTGAGAAAATAGACA R TTCTTCAGGAAAACACCTTCTTTGGACTCACA

Forage SNP @ 783 (unpadded: 734) P(SNP): 0.5604 P(VAR): 0.5603
       Freq3: 0.75 | B3: 41 B2: 33 | A: 0 C: 1 G: 0 T: 3
       TTTTGGGATAACCTGGATCCATAGATCGTTTA Y ATTCATCATACCTCCAGTATTTGTTAGCAACA
```

## Publication

Per Unneberg, Michael Strömberg, Fredrik Sterky, SNP discovery using advanced algorithms and neural networks, Bioinformatics, Volume 21, Issue 10, , Pages 2528–2530, https://doi.org/10.1093/bioinformatics/bti354
