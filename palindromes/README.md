# The analysis of metrics of palindromes accompanying the paper The Complete Sequence and Comparative Analysis of Ape Sex Chromosomes


We analyzed GC content of palindrome arms and a spacer, their respective lengths, and the identity of corresponding arms in each palindrome. 

### The analysis of correlations

The information in the main text, as well as Table S15 and plots from the Figure S7 were generated using the following R Sweave script:

```bash
palindrome_metrics.Rnw
```

If you plan to load the tables directly, and explore the data, please skip the first code chunk (the default behavior). This chunk merges the following datasets into a single table for X and Y chromosomes, respectively: 

```bash
palindromes.armlengths.txt
palindromes.spacerlengths.txt
palindroms.identity.txt
palindroms.geecee.txt
spacers.geecee.txt
palindromes.keyvaluepairs.txt 
spacers.keyvaluepairs.txt
```
These files should be available in two folders, one for each chromosome. If you are interested in how these files are generated, please continue reading below.

---

### Generating of the underlying data

In order to follow the results of this analysis, please start by creating a conda environment using file environment.yml

In the folder where you plan to run the analysis, please prepare the assemblies and the palindrome coordinates from palindrover. The files should be in the following format:

```bash
chrY.${species}.dip.20221111.fasta
${species}.chrY.palindromes.bed
chrX.${species}.dip.20221111.fasta
${species}.chrX.palindromes.bed 
```

using species names "mPanTro3" "mPanPan1" "mGorGor1" "mPonAbe1" "mPonPyg2" "mSymSyn1", for chimpanzee, bonobo, gorilla, Sumatran Orangutan, Bornean Orangutan, and siamang, respectively.

Now follow the instructions in the script:

```bash
./analyze_palindrome_metrics.sh
```
---
Contact <a href="mailto:cechova.biomonika@gmail.com">cechova.biomonika@gmail.com</a> with any questions. Please note that the code presented here is intended to support the aforementioned paper, not serve as a tool for the discovery and analysis of palindromes in new/additional species.
