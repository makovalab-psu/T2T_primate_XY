# This set of script annotates sequence classes on the X and Y chromosome of primates species and accompanies the paper The Complete Sequence and Comparative Analysis of Ape Sex Chromosomes

Per Skaletsky et al. 2003, ampliconic regions are long multi-copy regions with >50% sequence identity. However, these regions might sometimes be difficult to differentiate from satellite/heterochromatic regions, that are also present in multiple copies. For this reason, the first step in the identification of ampliconic regions, is hard-masking the X and Y chromosomes. Such hard-masking includes satellites, e.g. HSAT in human and subterminal satellites in bonobo, chimpanzee and gorilla. Additionally, features annotated as GAP and SAT are also hardmasked, as well as long interspersed elements (LINE).  

### Wide-masking of the chromosomes and mapping
This process requires repeatmasking primate chromosomes, spliting them into windows, and mapping them back to the chromosome assemblies in order to assess intrachromosomal similarity. It consists of three steps:

1. After this step only euchromatic portion of the Y chromosomes should remain
```
./repeatmask.sh 
```

2. Generates 5kb windows with 2kb steps
```
./make_windows.sh
```

2. Assigns the windows back to the assembly
```
./run_blastn.sh
```

Finally, windows of 5kb (using 2kb overlaps) along the remaining chromosomal sequences are mapped back to their respective assemblies. Long consecutive "hits" are then indicative of ampliconic regions. 


### Classification into the sequence classes

Please note the following considerations for the files used in the scripts below:

***merged.repeats.species.bed*** => merged coordinates of the continous stretches of repeat sequences. These long satellite stretches will be annotated as SATELLITE class. 

***ampliconic.bed*** => merged coordinates of the regions with high intrachromosomal similarity. These files are the output of the blastn window analysis described above. These regions will likely become part of AMPLICONIC class. 

***PARs.bed*** => coordinates of pseudoautozomal regions from Bob Harris. These regions will be annotated as PAR class.

XDEG_genes.bed and AMPL_genes.bed (for the Y chromosome only) => coordinates of X-degenerate and ampliconic genes from Karol Pal. These genes will help determine if regions will become ANCESTRAL/X-degenerate or AMPLICONIC.

***palindrover.merged.bed*** => coordinates of merged/flattened palindromes from Bob Harris. These regions will become part of AMPLICONIC class. 

Finally, the script tries to classify regions that would otherwise remained unclassified; for example, by "filling" in class information based on the adjacent regions (for example, a region nested within two ampliconic regions would also be considered ampliconic), or by merging the annotations of consecutive regions. The process of classification consists of the following steps:

1. run prepare_ampliconic_regions.sh (the script analyzes the results of the intrachromosomal similarity as revealed by blastn, and filters them) )
```
./prepare_ampliconic_regions.sh
```

2. The script generates both bed files for further downstream analysis, as well as files suitable for plotting with circos)
```
./create_sequence_classes_Y.sh
```

3. The script generates both bed files for further downstream analysis, as well as files suitable for plotting with circos)
```
./create_sequence_classes_X.sh
```

4. manually adjust XTR for both X and Y annotations


### Note on labeling

In the final publication, references to "X-Degenerate" (ancestral regions on
the Y chromosome) and "X-Ancestral" (ancestral regions on the X chromosome)
were changed to "Ancestral". Data and code herein uses "X-Degenerate" (or some
other form, e.g., "XDEG"). In the final sequence classes, released as
Additional File 2 with the manuscript, "XDEG" was renamed to "ANCESTRAL".
Similarly, "OTHER" was renamed to "ANCESTRAL" and "UNCLASSIFIED" on the X and Y
chromosomes, respectively.

---

Contact <a href="mailto:cechova.biomonika@gmail.com">cechova.biomonika@gmail.com</a> with any questions. Please note that the code presented here is intended to support the aforementioned paper, not serve as a tool for the sequence class annotation in new/additional species.
