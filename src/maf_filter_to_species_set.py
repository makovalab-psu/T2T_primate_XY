#!/usr/bin/env python3
"""
Read alignments in maf format and output those blocks that contain a specified
set of species, along with possibly other species; components for the other
species are discarded.
"""

from sys            import argv,stdin,stdout,stderr,exit
from math           import ceil
from maf_alignments import maf_alignments


programName    = "maf_filter_to_species_set"
programVersion = "0.1.0"


def usage(s=None):
	message = """

usage: cat maf | %s [options]
  --species=<list>       (cumulative) comma-separated list of species
  --pairwise             convert all surviving blocks to simple pairwise
                         blocks; note that no stitching is performed
  --discard:weeds        discard any components that go beyond the end of the
                         sequence ('off in the weeds')
                         (by default we treat these as errors and halt)
  --skip=<number>        skip the first so many maf blocks
  --head=<number>        limit the number of maf blocks read
  --progress=<number>    periodically report how many maf blocks we've
                         processed
  --version              report this program's version number

Read alignments in maf format and output those blocks that contain a specified
set of species. These blocks might have other species in the input, but
any components for the other species are discarded from the output.

The input is a maf file, as per the spec at
    http://genome.ucsc.edu/FAQ/FAQformat.html#format5. """ \
% programName

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	speciesOfInterest = None
	outputPairwise    = False
	discardWeeds      = False
	skipCount         = None
	headLimit         = None
	reportProgress    = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--species=")):
			if (speciesOfInterest == None):
				speciesOrder      = []
				speciesOfInterest = set()
			for species in argVal.split(","):
				species = species.strip()
				if (species not in speciesOfInterest): speciesOrder += [species]
				speciesOfInterest.add(species)
		elif (arg == "--pairwise"):
			outputPairwise = True
		elif (arg in ["--discard:weeds","--discard=weeds"]):
			discardWeeds = True
		elif (arg.startswith("--skip=")):
			skipCount = int_with_unit(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg in ["--version","--v","--V","-version","-v","-V"]):
			exit("%s, version %s" % (programName,programVersion))
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (speciesOfInterest == None):
		usage("a species list is required")

	# process the alignments

	needHeaderSpacer  = False
	nothingHasPrinted = True

	skipsRemaining = skipCount
	mafBlockNumber = mafBlocksSkipped = 0
	mafBlocksKept = mafBlocksModified = mafBlocksDiscarded = 0
	for a in maf_alignments(stdin,yieldHeader=True,saveLines=True,discardWeeds=discardWeeds):
		if (type(a) == str):  # this is a maf header
			print(a)
			needHeaderSpacer = True
			continue

		if (needHeaderSpacer):
			print()
			needHeaderSpacer = False

		if (skipsRemaining != None):
			mafBlocksSkipped += 1
			if (mafBlocksSkipped <= skipsRemaining): continue
			print("first %s maf blocks skipped" % "{:,}".format(skipCount),file=stderr)
			skipsRemaining = None

		mafBlockNumber += 1
		if (headLimit != None) and (mafBlockNumber > headLimit):
			print("limit of %s input maf blocks reached" % "{:,}".format(headLimit),file=stderr)
			break
		if (reportProgress != None) and (mafBlockNumber % reportProgress == 0):
			print("progress: maf block %s (line %s), %s kept as is, %s modified, %s discarded" \
			    % ("{:,}".format(mafBlockNumber),
			       "{:,}".format(a.lineNum),
			       "{:,}".format(mafBlocksKept),
			       "{:,}".format(mafBlocksModified),
			       "{:,}".format(mafBlocksDiscarded)),
			      file=stderr)

		refs = set([c.ref for c in a.block])
		if (not speciesOfInterest.issubset(refs)):
			mafBlocksDiscarded += 1
			continue

		if (refs == speciesOfInterest):
			mafBlocksKept += 1
		else:
			mafBlocksModified += 1

		if (outputPairwise):
			refToComponents = {name:[] for name in speciesOfInterest}
			for c in a.block:
				if (c.ref in speciesOfInterest): refToComponents[c.ref] += [c]
			for (name1,name2) in all_pairs(speciesOrder,ordered=False):
				for c1 in refToComponents[name1]:
					for c2 in refToComponents[name2]:
						if (nothingHasPrinted): nothingHasPrinted = False
						else:                   print()
						print("a")
						print(c1.line)
						print(c2.line)
		else:
			if (nothingHasPrinted): nothingHasPrinted = False
			else:                   print()
			print("a")
			for c in a.block:
				if (c.ref in speciesOfInterest): print(c.line)

	if (reportProgress != None) and (headLimit == None):
		print("progress: %s maf blocks total" % ("{:,}".format(mafBlockNumber)),file=stderr)

	if (mafBlocksKept + mafBlocksModified + mafBlocksDiscarded == 0):
		print("no maf blocks kept, none discarded",file=stderr)
	else:
		print("%s maf blocks kept as is (%.2f%%), %s modified (%.2f%%), %s discarded" \
		    % ("{:,}".format(mafBlocksKept),
		       100.0*mafBlocksKept/(mafBlocksKept+mafBlocksModified+mafBlocksDiscarded),
		       "{:,}".format(mafBlocksModified),
		       100.0*mafBlocksModified/(mafBlocksKept+mafBlocksModified+mafBlocksDiscarded),
		       "{:,}".format(mafBlocksDiscarded)),
		      file=stderr)


# all_pairs--
#	yield all pairs from a set
#
# if ordered is false, pairs (a,b) and (b,a) are considered the same, and only
# one is generated

def all_pairs(items,ordered=False):
	for ix in range(len(items)-1):
		for iy in range(ix+1,len(items)):
			yield (items[ix],items[iy])
			if (ordered): yield (items[iy],items[ix])


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
