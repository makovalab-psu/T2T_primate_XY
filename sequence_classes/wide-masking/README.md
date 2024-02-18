## This page describes how to repeatmask primate Y chromosomes, split them into windows, and map back to the Y chromosome assemblies in order to assess intrachromosomal similarity

Run the scripts in the following order:
1) repeatmask.sh (after this step only euchromatic portion of the Y chromosomes should remain)
2) make_windows.sh (generates 5kb windows with 2kb steps)
3) run_blastn.sh (assigns the windows back to the assembly)
