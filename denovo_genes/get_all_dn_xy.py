from Bio import AlignIO

####bonobo de novo gene regions

idx = AlignIO.MafIO.MafIndex("panpanY.mafindex", "mPanPan1.dip.20221111.chrY_onto_X.maf", "mPanPan1.chrY")
multiple_alignment1 = idx.get_spliced(
    [28075448,28076422,28077214,28077861],
    [28075528,28076556,28077296,28078641
],
    1
)

idx = AlignIO.MafIO.MafIndex("pantroY.mafindex", "mPanTro3.dip.20221111.chrY_onto_X.maf", "mPanTro3.chrY")
multiple_alignment2 = idx.get_spliced(
    [29340921,29341894,29342686,29343339],
    [29341001,29342030,29342769,29344121],
    1
)

idx = AlignIO.MafIO.MafIndex("humanY.mafindex", "chrY_hg002_v2.7_onto_chm13.chrX.maf", "hg002_v2.chrY")
multiple_alignment3 = idx.get_spliced(
    [10515336,10517427,10516688,10518456],
    [10516118,10517563,10516771,10518536],
    -1
)

bonobo_alignments=[multiple_alignment1,multiple_alignment2, multiple_alignment3]
AlignIO.write(bonobo_alignments, "panpan_new_dn_xy_aln_all.fa", "fasta")


###siamang de novo gene regions

idx = AlignIO.MafIO.MafIndex("symsynY.mafindex", "mSymSyn1.dip.20221111.chrY_onto_X.maf", "mSymSyn1.chrY")
multiple_alignment4 = idx.get_spliced(
    [11675362,11646093,11640574,11594759,11592382,11591836],
    [11675858,11646226,11640885,11594899,11592503,11592280],
    -1
)

idx = AlignIO.MafIO.MafIndex("gorgorX.mafindex", "mGorGor1.dip.20221111.chrY_onto_X.maf", "mGorGor1.chrX")
multiple_alignment5 = idx.get_spliced(
    [17622742,17688219,17728182,17732272,17732542],
    [17623104,17688381,17728327,17732399,17733019],
    1
)

idx = AlignIO.MafIO.MafIndex("ponpygX.mafindex", "mPonPyg2.dip.20221111.chrY_onto_X.maf", "mPonPyg2.chrX")
multiple_alignment6 = idx.get_spliced(
    [7535980,7639756,7643882,7644151],
    [7536335,7639901,7644009,7644614],
    1
)

idx = AlignIO.MafIO.MafIndex("ponabeX.mafindex", "mPonAbe1.dip.20221111.chrY_onto_X.maf", "mPonAbe1.chrX")
multiple_alignment7 = idx.get_spliced(
    [7128374,7233595,7237722,7237991],
    [7128729,7233740,7237849,7238462],
    1
)

siamang_alignments=[multiple_alignment4,multiple_alignment5, multiple_alignment6, multiple_alignment7]
AlignIO.write(siamang_alignments, "symsyn_new_dn_xy_aln_all.fa", "fasta")


