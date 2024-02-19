#!/usr/bin/env python3
"""
Given pairwise self-alignments in lastz's general format, identify palindromes
"""

from sys                import argv,stdin,stdout,stderr,exit
from math               import ceil
from collections        import defaultdict
from copy               import copy
from gzip               import open as gzip_open
from alignment_table    import AlignmentTable
from cigar              import split_cigar,cigar_to_string,cigar_identity,cigar_nmatch
from transitive_closure import TransitiveClosure


programName    = "palindrover"
programVersion = "0.1.5"


def usage(s=None):
	message = """
usage: cat <alignment_file> | %s [options]
  --[min]identity=<value>    (ID=) ignore alignments with identity below
                             <value>; 0<=value<=1, and a percent sign may be
                             used
                             (default is 98%%)
  --[min]length=<number>     (L=) minimum length of one arm of a palindrome
                             (default is 100K)
  --[max]spacer=<number>     (S=) maximum length of the spacer between
                             palindrome arms
                             (default is 500K)
  --[min]identity:reject=<value> minimum identity for a rejected
                             palindrome for it to be reported to the rejects
                             file
                             (default is whatever is used for palindromes)
  --[min]length:reject=<number> minimum length of one arm of a rejected
                             palindrome for it to be reported to the rejects
                             file
                             (default is 1)
  --blacklist=<file>         (cumulative) file containing 'undesirable'
                             intervals; candidates that are blacklisted are
                             are not reported, but are *not* considered to be
                             rejects
                             (see file format below)
  --blacklist[:<threshold>]=<file>  same as --blacklist=<file> but a candidate
                             is blacklisted only if the portion of its arm
                             overlapping 'undesirable' intervals exceeds the
                             threshold; <threshold> is in the range 0..1 and
                             can be expressed as a number, percentage, or
                             fraction
  --report:blacklisted       report blacklisted alignments (to stderr)
  --column:blacklisted       add a "blacklisted%%" column to the output
  --column:palname           add a "palName" column to the output; note that we
                             assign names like "Q1" and "Q2" so that they will
                             not be confused with names like "P1" and "P2" in
                             the literature.
  --group:overlaps           group palindromes that overlap; this *only*
                             affects the names in the palName column, as
                             palindromes in the same overlap group will be
                             given names like "Q1.1" and "Q1.2"
  --alias:<alias>=<name>     the input can use <alias> as an alias for column
                             name <name>
  --output=<filename>        write palindromes to a file
                             (by default, we write palindromes to stdout)
  --output:bed=<filename>    additionally, write palindrome intervals to a
                             three- or four-column bed file
  --rejects=<filename>       write rejected candidates to a file
                             (by default, we don't report rejected candidates)
  --blacklisteds=<filename>  write blacklisted candidates to a file
  --head=<number>            limit the number of alignment records
  --progress=<number>        periodically report how many alignment records
                             we've  processed
  --version                  report this program's version number

Read pairwise self-alignments from lastz's general format and identify
palindromes.

An underlying premise is that self-alignment is a symmetric process. Thus we
assume any alignment on the reverse strand has a twin mirrored on the other
side of the main diagonal. We do not check this and simply operate as though it
is true; so we only look for alignments above the main diagonal.

Any alignments of different sequences are ignored, as are any alignments on the
forward (+) strand, or any portion of an alignment that is below the main
diagonal.

Typical input is shown below but other columns may be included, and columns can
appear in any order. Cigarx is required; specifically, we need cigar in the
format that distinguishes between matches and mismatches (operators in
{=,X,I,D}).

Such input can be created by lastz with a command like this:
  lastz mule.fa mule.fa --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%%,cigarx
However, any aligner can be used, so long as it produces the indicated fields.
The specifics of those fields can be found at
  https://lastz.github.io/lastz/#fmt_general

Typical input:
  #name1 zstart1 end1   name2 strand2 zstart2+ end2+   id%%  cigarx
  chrY   108173  108756 chrY  +       107993   108576  89.0%% ...
  chrY   127815  127994 chrY  -       3472294  3472475 85.9%% ...
  chrY   472065  472284 chrY  +       472213   472432  82.6%% ...
  chrY   592387  593300 chrY  +       592697   593620  91.4%% ...
  chrY   606857  606925 chrY  -       606857   606925  97.1%% ...
  chrY   867024  867083 chrY  -       866643   866702  98.3%% ...
  chrY   917166  917759 chrY  -       918863   919458  88.4%% ...
  chrY   959027  959803 chrY  +       957987   958763  84.8%% ...
   ...

Output is the same format as the input, but additional columns palArmLen and
palSpacer are added, as are, optionally, blacklisted%% and palName columns.

If a <threshold> is expressed for any blacklist file the same threshold is then
applied for all files. If no threshold is specified then any intersection with
a blacklist interval disqualifies an alignment.

Blacklist files have lines of the form
    <chrom> <start> <end> [<extras>]
where <start> and <end> are origin-zero half-open. Anything beyond the third
column is ignored.

If --output:bed is specified, a bed file is written, per the spec at
    http://genome.ucsc.edu/FAQ/FAQformat.html#format1
Only the chrom, start, and end fields are written, and optionally name. A line
is written for each palindrome arm, with the suffixes "A" and "B" added to the
name to distinguish them.""" \
% programName

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


requiredColumns    = ("name1","zstart1","end1",
                      "name2","zstart2+","end2+",
                      "strand",
                      "id","cigarx")
nonRequiredColumns = ("nmatch",)                 # the comma makes this a tuple
columnAliases      = {"strand2" : "strand",
                      "s"       : "strand",
                      "s2"      : "strand",
                      "id%"     : "id"}


def main():
	global debug

	# parse the command line

	minIdentity          = 0.98
	minLength            = 100*1000    # length of one arm
	maxSpacer            = 500*1000    # distance between the arms
	minRejectIdentity    = None
	minRejectLength      = None
	blacklistThreshold   = None
	blacklistFilenames   = []
	reportBlacklisted    = False
	createBlacklistedColumn = False
	createPalNameColumn  = False
	groupOverlaps        = False
	aliasToName          = {}
	outputFilename       = None
	bedFilename          = None
	rejectsFilename      = None
	blacklistedsFilename = None
	headLimit            = None
	reportProgress       = None
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if ((arg.startswith("--identity="))
		 or (arg.startswith("--minidentity="))
		 or (arg.startswith("--id="))
		 or (arg.startswith("--minid="))
		 or (arg.startswith("ID="))):
			minIdentity = parse_probability(argVal)
		elif ((arg.startswith("--length="))
		 or (arg.startswith("--minlength="))
		 or (arg.startswith("L="))):
			minLength = int_with_unit(argVal)
		elif ((arg.startswith("--spacer="))
		 or (arg.startswith("--maxspacer="))
		 or (arg.startswith("S="))):
			maxSpacer = int_with_unit(argVal)
		elif ((arg.startswith("--identity:reject="))
		 or (arg.startswith("--minidentity:reject="))
		 or (arg.startswith("--id:reject="))
		 or (arg.startswith("--minid:reject="))):
			minRejectIdentity = parse_probability(argVal)
		elif ((arg.startswith("--length:reject="))
		 or (arg.startswith("--minlength:reject="))):
			minRejectLength = int_with_unit(argVal)
		elif (arg.startswith("--blacklist=")):
			blacklistFilenames += [argVal]
		elif (arg.startswith("--blacklist:")):
			argVal = arg.split(":",1)[1]
			if ("=" not in argVal):  # --blacklist:<threshold>
				argVal = arg.split(":",1)[1]
				t = parse_probability(argVal)
				if (blacklistThreshold == None):
					blacklistThreshold = t
				elif (t != blacklistThreshold):
					usage("blacklist threshold %0.5f (in %s) disagrees with %0.5f specified earlier" \
					    % (t,arg,blacklistThreshold))
			else:                    # --blacklist:<threshold>=<file>
				(t,filename) = argVal.split("=",1)
				t = parse_probability(t)
				if (blacklistThreshold == None):
					blacklistThreshold = t
				elif (t != blacklistThreshold):
					usage("blacklist threshold %0.5f (in %s) disagrees with %0.5f specified earlier" \
					    % (t,arg,blacklistThreshold))
				blacklistFilenames += [filename]
		elif (arg in ["--report:blacklisted","--report:blacklist"]):
			reportBlacklisted = True
		elif (arg in ["--column:blacklisted","--column:blacklisted%"]):
			createBlacklistedColumn = True
		elif (arg in ["--column:name","--column:palname","--column:palName"]):
			createPalNameColumn = True
		elif (arg == "--group:overlaps"):
			groupOverlaps = True
		elif (arg.startswith("--alias:")):
			argVal = arg.split(":",1)[1]
			for argField in argVal.split(","):
				(alias,name) = argField.split("=",1)
				assert (alias not in aliasToName)
				aliasToName[alias] = name
		elif (arg.startswith("--output=")) or (arg.startswith("--out=")):
			outputFilename = argVal
		elif (arg.startswith("--output:bed=")) or (arg.startswith("--out:bed=")):
			bedFilename = argVal
		elif (arg.startswith("--rejects=")):
			rejectsFilename = argVal
		elif (arg.startswith("--blacklisteds=")) or (arg.startswith("--blacklisted=")):
			blacklistedsFilename = argVal
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

	if (blacklistThreshold != None) and (blacklistFilenames == []):
		usage("a blacklist threshold was specified, but no --blacklist files were provided")

	if (reportBlacklisted) and (blacklistFilenames == []):
		usage("--report:blacklisted was specified, but no --blacklist files were provided")

	if (rejectsFilename == None):
		if (minRejectIdentity != None):
			usage("--[min]identity:reject was specified, but no --rejects file to write to")
		if (minRejectLength != None):
			usage("--[min]length:reject was specified, but no --rejects file to write to")
	else: # if (rejectsFilename != None):
		if (minRejectIdentity == None): minRejectIdentity = minIdentity
		if (minRejectLength   == None): minRejectLength = 0

	if ("blacklisted%" in debug):   # (for backwards compatibility)
		createBlacklistedColumn = True

	if (groupOverlaps) and (not createPalNameColumn):
		groupOverlaps = False
		print("--group:overlaps is ignored, since --column:palname was not specified",file=stderr)

    # read the blacklist files, and collapse overlaps

	blacklistingEnabled = (blacklistFilenames != [])

	chromToBlacklist = {}

	for filename in blacklistFilenames:
		if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
			f = gzip_open(filename,"rt")
		else:
			f = open(filename,"rt")
		for (lineNumber,chrom,start,end) in read_intervals(f,filename):
			if (chrom not in chromToBlacklist): chromToBlacklist[chrom] = []
			chromToBlacklist[chrom] += [(start,end)]
		f.close()

	for chrom in chromToBlacklist:
		chromToBlacklist[chrom] = merge_overlapping_intervals(chromToBlacklist[chrom])

	# open the alignment table, as an iterator

	for alias in aliasToName:
		columnAliases[alias] = aliasToName[alias]

	t = AlignmentTable.from_file(stdin,
	                             requiredColumns=requiredColumns,
	                             nonRequiredColumns=nonRequiredColumns,
	                             columnAliases=columnAliases)
	header = t.header[1:].split()
	nameToColumn = dict(t.columnNames)

	# rearrange output columns so that if cigarx was the final column, it
	# remains as such

	newColumns = ["palArmLen","palSpacer"]
	if (createPalNameColumn):
		newColumns = ["palName"] + newColumns
	if (createBlacklistedColumn):
		newColumns = ["blackListed%"] + newColumns

	if (header[-1] == "cigarx"):
		header = header[:-1] + newColumns + ["cigarx"]
		for i in range(1,len(newColumns)+2):
			nameToColumn[header[-i]] = len(header)-i
	else:
		header += newColumns
		for i in range(1,len(newColumns)+1):
			nameToColumn[header[-i]] = len(header)-i

	# open/create files

	outF = stdout
	if (outputFilename != None):
		if (outputFilename.endswith(".gz")) or (outputFilename.endswith(".gzip")):
			outF = gzip_open(outputFilename,"wt")
		else:
			outF = open(outputFilename,"wt")

	rejectsF = rejects = None
	if (rejectsFilename != None):
		rejects = []
		if (rejectsFilename.endswith(".gz")) or (rejectsFilename.endswith(".gzip")):
			rejectsF = gzip_open(rejectsFilename,"wt")
		else:
			rejectsF = open(rejectsFilename,"wt")

	blacklistedsF = blacklisteds = None
	if (blacklistedsFilename != None):
		blacklisteds = []
		if (blacklistedsFilename.endswith(".gz")) or (blacklistedsFilename.endswith(".gzip")):
			blacklistedsF = gzip_open(blacklistedsFilename,"wt")
		else:
			blacklistedsF = open(blacklistedsFilename,"wt")

	bedF = None
	if (bedFilename != None):
		if (bedFilename.endswith(".gz")) or (bedFilename.endswith(".gzip")):
			bedF = gzip_open(bedFilename,"wt")
		else:
			bedF = open(bedFilename,"wt")

	# process the alignments, collecting palindromes (and rejects)

	minSearchIdentity = minIdentity if (minRejectIdentity == None) else min(minIdentity,minRejectIdentity)

	palindromes = []
	numBlacklisted = 0
	alignmentNumber = 0
	for a in t:
		alignmentNumber += 1

		if (headLimit != None) and (alignmentNumber > headLimit):
			print("limit of %s alignments reached" % ("{:,}".format(headLimit)),file=stderr)
			break
		if (reportProgress != None):
			if (alignmentNumber == 1) or (alignmentNumber % reportProgress == 0):
				print("processing alignment %s" % "{:,}".format(alignmentNumber),file=stderr)

		if (a.name1 != a.name2): continue
		if (a.strand != "-"):    continue

		aPal = palindrome(a,minSearchIdentity,minLength,maxSpacer,reportRejects=rejectsF!=None)
		if (aPal == None): continue

		if (createBlacklistedColumn):
			bPct = fraction_blacklisted(aPal,chromToBlacklist)
			aPal.fractionBlacklisted = max(bPct.leftFractionCovered,bPct.rightFractionCovered)

		if (aPal.isPalindrome):
			if (aPal.armLength < minLength): aPal.isPalindrome = False
		if (aPal.isPalindrome) and (minSearchIdentity < minIdentity):
			if (aPal.idFloat < minIdentity): aPal.isPalindrome = False

		if (not aPal.isPalindrome):
			if (aPal.armLength < minRejectLength): continue
		if (not aPal.isPalindrome) and (minSearchIdentity < minRejectIdentity):
			if (aPal.idFloat < minRejectIdentity): continue

		bInfo = is_blacklisted(aPal,chromToBlacklist,blacklistThreshold)
		if (bInfo != None):
			if (blacklisteds != None): blacklisteds += [aPal]
			numBlacklisted += 1
			if (not reportBlacklisted):   continue
			if (aPal.start1 > aPal.end2): continue    # (don't report if below main diagonal)
			if (blacklistThreshold == None):
				(s,e) = bInfo.interval
				print("(at line %s) %s %d..%d %d..%d discarded; blacklist %d..%d" \
				    % ("{:,}".format(aPal.lineNumber),aPal.name1,aPal.start1,aPal.end1,aPal.start2,aPal.end2,s,e),
				      file=stderr)
			else: # if (blacklistThreshold != None):
				print("(at line %s) %s %d..%d %d..%d discarded; blacklist %.1f%% (%s part)" \
				    % ("{:,}".format(aPal.lineNumber),aPal.name1,aPal.start1,aPal.end1,aPal.start2,aPal.end2,100*bInfo.fractionCovered,bInfo.side),
				      file=stderr)
			continue

		if (aPal.isPalindrome):
			palindromes += [aPal]
		elif (rejects != None):
			rejects += [aPal]

	# deal with "nested" overlaps in the palindrome set

	if (groupOverlaps):
		clusterer = cluster_by_overlaps(palindromes)
		palIxToGroupMembers = clusterer.derive_clusters()
		for palIx in palIxToGroupMembers:
			members = list(palIxToGroupMembers[palIx])
			members.sort()
			minPalIx = members[0]
			for (memberId,memberPalIx) in enumerate(members):
				aPal = palindromes[memberPalIx]
				aPal.groupId  = minPalIx
				aPal.memberId = memberId+1

	# report palindromes (and rejects)

	print("#%s" % "\t".join(header),file=outF)
	nameToPalindromNumber = defaultdict(int)
	groupIdToName = {}
	for aPal in palindromes:
		if (createPalNameColumn):
			if (not groupOverlaps) or (not hasattr(aPal,"groupId")):
				nameToPalindromNumber[aPal.name1] += 1
				aPal.palName = "%s.Q%d" % (aPal.name1,nameToPalindromNumber[aPal.name1])
			elif (aPal.memberId == 1):
				nameToPalindromNumber[aPal.name1] += 1
				groupIdToName[aPal.groupId] = "%s.Q%d" % (aPal.name1,nameToPalindromNumber[aPal.name1])
				aPal.palName = "%s.1" % groupIdToName[aPal.groupId]
			else:
				aPal.palName = "%s.%d" % (groupIdToName[aPal.groupId],aPal.memberId)
		print(create_alignment_line(aPal,nameToColumn),file=outF)
		if (bedF != None):
			print(create_bed_lines(aPal),file=bedF)

	if (rejects != None):
		print("#%s" % "\t".join(header),file=rejectsF)
		nameToPalindromNumber = defaultdict(int)
		for aPal in rejects:
			if (createPalNameColumn):
				nameToPalindromNumber[aPal.name1] += 1
				aPal.palName = "%s.R%d" % (aPal.name1,nameToPalindromNumber[aPal.name1])
			print(create_alignment_line(aPal,nameToColumn),file=rejectsF)

	if (blacklisteds != None):
		print("#%s" % "\t".join(header),file=blacklistedsF)
		nameToPalindromNumber = defaultdict(int)
		for aPal in blacklisteds:
			if (createPalNameColumn):
				nameToPalindromNumber[aPal.name1] += 1
				aPal.palName = "%s.B%d" % (aPal.name1,nameToPalindromNumber[aPal.name1])
			print(create_alignment_line(aPal,nameToColumn),file=blacklistedsF)

	if (outputFilename       != None): outF.close()
	if (rejectsFilename      != None): rejectsF.close()
	if (blacklistedsFilename != None): blacklistedsF.close()
	if (bedFilename          != None): bedF.close()

	if (blacklistingEnabled):
		if (numBlacklisted == 0):   print("no palindromes were blacklisted",file=stderr)
		elif (numBlacklisted == 1): print("1 palindromes was blacklisted",file=stderr)
		else:                       print("%s palindromes were blacklisted" % "{:,}".format(numBlacklisted),file=stderr)


# is_blacklisted--
#	Determine whether the portion of an alignment that intersects a blacklisted
#	interval exceeds a specified threshold.
#
#	Returns None if the alignment is not blacklisted. Otherwise it returns an
#	object with information about why it is blacklisted.

class BlacklistInfo: pass

def is_blacklisted(a,chromToBlacklist,blacklistThreshold):
	assert (a.name1 == a.name2)
	assert (a.strand == "-")

	if (a.name1 not in chromToBlacklist): return None

	intervals = chromToBlacklist[a.name1]

	# if there is no threshold, we just report the first overlapping interval
	# we find

	if (blacklistThreshold == None):
		blacklistInterval = find_overlap(intervals,a.start1,a.end1)
		if (blacklistInterval == None):
			blacklistInterval = find_overlap(intervals,a.start2,a.end2)

		if (blacklistInterval == None): return None

		bInfo = BlacklistInfo()
		bInfo.interval = blacklistInterval
		return bInfo

	# otherwise, we compute the portion overlapping, separately for each 'side'
	# of the alignment; if either side exceeds the threshold we report it

	blacklistCount = count_overlap_bases(intervals,a.start1,a.end1)
	if (blacklistCount >= blacklistThreshold * (a.end1-a.start1)):
		bInfo = BlacklistInfo()
		bInfo.side            = "left"
		bInfo.fractionCovered = float(blacklistCount) / (a.end1-a.start1)
		return bInfo

	blacklistCount = count_overlap_bases(intervals,a.start2,a.end2)
	if (blacklistCount >= blacklistThreshold * (a.end2-a.start2)):
		bInfo = BlacklistInfo()
		bInfo.side            = "right"
		bInfo.fractionCovered = float(blacklistCount) / (a.end2-a.start2)
		return bInfo

	return None


# fraction_blacklisted--
#	Compute the portion of a palindrome that is blacklisted.

def fraction_blacklisted(a,chromToBlacklist):
	assert (a.name1 == a.name2)
	assert (a.strand == "-")

	if (a.name1 not in chromToBlacklist): return 0.0

	intervals = chromToBlacklist[a.name1]

	# otherwise, we compute the portion overlapping, separately for each 'side'
	# of the alignment; if either side exceeds the threshold we report it

	bInfo = BlacklistInfo()

	blacklistCount = count_overlap_bases(intervals,a.start1,a.end1)
	bInfo.leftFractionCovered = float(blacklistCount) / (a.end1-a.start1)

	blacklistCount = count_overlap_bases(intervals,a.start2,a.end2)
	bInfo.rightFractionCovered = float(blacklistCount) / (a.end2-a.start2)

	return bInfo


# palindrome--
#	Determine whether an alignment is the left arm of a palindrome. If it is,
#	return a modified copy of the alignment, represent the palindrome. Otherwise
#	return None
#
# note that minRejectIdentity and minRejectLength are *not* handled by this
# function

def palindrome(a,minIdentity,minLength,maxSpacer,reportRejects=False):
	# comments herein relate to a dotplot with position 1 increasing to the
	# east and position 2 increasing to the north; also note that start2,end2
	# is counted along forward strand, but alignment is on reverse

	assert (a.name1 == a.name2)
	assert (a.strand == "-")

	identity = float(a.id[:-1])/100.0 if (a.id.endswith("%")) else float(a.id)
	if (identity < minIdentity): return None

	if (a.start1 > a.end2): return None                      # nw end is below main diagonal

	isReject = False

	aMod = None
	if (a.end1 > a.start2):                                  # se end is below main diagonal
		armLength = (a.end2-a.start1) / 2.0
		if (armLength < minLength) and (not reportRejects):  # alignment is too short
			return None
		aMod = copy(a)
		truncate_cigar_at_diagonal(aMod)                     # truncate, then recompute
		twiceDistToMainDiagonal = aMod.start2 - aMod.end1
		armLength = aMod.end1 - aMod.start1
		spacerLength = twiceDistToMainDiagonal
	else:                                                    # se end is above main diagonal
		twiceDistToMainDiagonal = a.start2 - a.end1
		armLength = a.end1 - a.start1
		spacerLength = twiceDistToMainDiagonal

	if (armLength < minLength):    isReject = True           # alignment is too short
	if (spacerLength > maxSpacer): isReject = True           # spacer is too long

	if (isReject) and (not reportRejects):
		return None

	if (aMod == None): aMod = copy(a)
	aMod.isPalindrome = not isReject
	aMod.idFloat      = identity
	aMod.armLength    = armLength
	aMod.spacerLength = spacerLength
	return aMod


# truncate_cigar_at_diagonal--
#	Follow an alignment's cigarx string and truncate it when it is about to
#	cross the main diagonal. The alignment's cigar string is modified, as are
#	start and ends
#
# Note that we may truncate at a point above the main diagonal. This can happen
# if, for instance, we reach a point 1 above the diagonal and the next
# operation is a match or mismatch (would take us past the diagaonal).
# Morevover, we will not include terminal indels or mismatches in the truncated
# list, and this may also leave us above the diagonal.

def truncate_cigar_at_diagonal(a):
	assert (a.strand == "-")
	assert (a.start1 <= a.end2)   # require nw end above (or on) main diagonal
	assert (a.end1 > a.start2)    # require se end below main diagonal

	cigarInfo = split_cigar(a.cigarx)

	# trace the path until we reach the main diagonal

	(tPos,qPos) = (a.start1,a.end2)
	truncatedOperations = []
	for (rpt,op) in cigarInfo.operations:
		# $$$ these two conditions can obviously be combined, but stating them
		#     separately seems clearer
		if (tPos >= qPos):               # we are at or below main diagonal,
			break                        # .. so we're done
		if (tPos == qPos-1):             # we are one above main diagonal,
			break                        # .. so we're also done

		if (op in ["M","X","="]):
			if (tPos+rpt > qPos-rpt):    # full rpt would take us below main diagonal
				rpt = (qPos-tPos) // 2
			tPos += rpt
			qPos -= rpt
		elif (op == "I"):
			if (tPos > qPos-rpt):        # full rpt would take us below main diagonal
				rpt = qPos-tPos
			qPos -= rpt
		elif (op == "D"):
			if (tPos+rpt > qPos):        # full rpt would take us below main diagonal
				rpt = qPos-tPos
			tPos += rpt
		else:
			raise ValueError

		if (rpt > 0): truncatedOperations += [(rpt,op)]

	# remove any series of indels or mismatches from the end of our truncated
	# path

	goodLen = len(truncatedOperations)
	while (goodLen > 1):
		(rpt,op) = truncatedOperations[goodLen-1]
		if (op not in ("I","D","X")): break
		if (op == "X"):
			tPos -= rpt
			qPos += rpt
		elif (op == "I"):
			qPos += rpt
		else: # if (op == "D"):
			tPos -= rpt
		goodLen -= 1

	truncatedOperations = truncatedOperations[:goodLen]

	# modify the alignment object

	a.cigarx = cigar_to_string(truncatedOperations)
	a.end1   = tPos
	a.start2 = qPos

# create_alignment_line--
#	Create a line, consistent with the input lines, for an alignment that has
#	been modified

def create_alignment_line(a,nameToColumn,fieldJoiner="\t"):
	# (we safely assume that columns and names are distinct)
	columnToName = {col:name for (name,col) in nameToColumn.items()}
	maxCol = max(columnToName.keys())

	cigarInfo = None

	line = []
	for col in range(maxCol+1):
		if (col not in columnToName):   # shouldn't happen
			line += ["?"]
			continue
		name = columnToName[col]
		if (name == "zstart1"):
			line += [str(a.start1)]
		elif (name == "end1"):
			line += [str(a.end1)]
		elif (name == "zstart2+"):
			line += [str(a.start2)]
		elif (name == "end2+"):
			line += [str(a.end2)]
		elif (name in ["id%","id"]):
			if (cigarInfo == None):
				cigarInfo = split_cigar(a.cigarx)
			line += ["%.3f%%" % (100*cigar_identity(cigarInfo))]
		elif (name == "nmatch"):
			if (cigarInfo == None):
				cigarInfo = split_cigar(a.cigarx)
			line += [str(cigar_nmatch(cigarInfo))]
		elif (name == "palArmLen"):
			line += [str(a.armLength)]
		elif (name == "palSpacer"):
			line += [str(a.spacerLength)]
		elif (name == "palName"):
			line += [a.palName]
		elif (name == "blackListed%"):
			line += ["%.1f%%" % (100*a.fractionBlacklisted)]
		elif (name in requiredColumns):
			line += [str(getattr(a,name))]
		# elif (name == "score"):
		#	to do score we'd need text1,text2, and to recreate those we'd need
		#   the underlying sequences
		else:
			line += ["?"]

	return fieldJoiner.join(line)


# create_bed_lines--
#	Create two lines, in bed format, describing the two arms of the palindrome.

def create_bed_lines(a):
	lines = []

	line  = [a.name1,str(a.start1),str(a.end1)]
	if (hasattr(a,"palName")): line += [a.palName+"A"]
	lines += ["\t".join(line)]

	line  = [a.name2,str(a.start2),str(a.end2)]
	if (hasattr(a,"palName")): line += [a.palName+"B"]
	lines += ["\t".join(line)]

	return "\n".join(lines)


# cluster_by_overlaps--
#	Gather and overlapping palindromes into groups. In the resulting
#	clustering object, item ids are indexes into the input list.

def cluster_by_overlaps(palindromes):
	# collect intervals by chromosome

	chromToPalindromes = defaultdict(list)
	for (palIx,aPal) in enumerate(palindromes):
		chromToPalindromes[aPal.name1] += [(aPal.start1,aPal.end1,palIx)]
		chromToPalindromes[aPal.name1] += [(aPal.start2,aPal.end2,palIx)]

	# find overlapping intervals and 'connect' them

	clusterer = TransitiveClosure()
	for chrom in chromToPalindromes:
		intervals = chromToPalindromes[chrom]
		intervals.sort()

		# scan the sorted intervals; ix points to the an interval that does not
		# overlap any earlier interval; ix2 points to an interval that overlaps
		# ix, scans forward until a non-overlapping interval is encountered, at
		# which point it points to the latest overlapping interval; as we scan
		# for ix2, each overlapping interval is 'connected' to ix, and because
		# the clusterer implements transitive closure, all the overlapping
		# intervals will end up in a cluster with ix

		ix = 0
		while (ix < len(intervals)):
			(start,end,palIx) = intervals[ix]
			ix2 = ix
			while (ix2+1 < len(intervals)):
				(start2,end2,palIx2) = intervals[ix2+1]
				if (start2 >= end): break   # (no overlap)
				end = max(end,end2)
				clusterer.connect(palIx,palIx2)
				if ("cluster" in debug):
					print("connect %d-%d [%d]-[%d] (%d..%d)-(%d..%d)" \
					    % (ix,ix2+1,palIx,palIx2,start,end,start2,end2),file=stderr)
				ix2 += 1
			ix = ix2+1

	return clusterer


# read_intervals--
#	Read a three-column file (ignoring extra columns) and yield each interval.

def read_intervals(f,filename):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		       "not enough fields at line %d in %s (%d, expected at least %d)" \
		     % (lineNumber,filename,len(fields),3)

		try:
			ref   = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (end < start): raise ValueError
		except ValueError:
			assert (False), "bad line (%d in %s): %s" % (lineNumber,filename,line)

		yield (lineNumber,ref,start,end)


# merge_overlapping_intervals--
#	collapse overlapping intervals and return intervals covering the same set
#	of positions, but with no overlaps; the resulting list is ordered by
#	increasing position.
#
# side effect: the incoming list is modified (it is sorted)

def merge_overlapping_intervals(intervals):
	if (intervals == []): return []

	intervals.sort()

	overlaps = []
	start = end = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s < end):
			end = max(end,e)
		else:
			overlaps += [(start,end)]
			(start,end) = (s,e)

	if (start != None):
		overlaps += [(start,end)]

	return overlaps


# find_overlap--
#	Determine whether the specified interval overlaps any interval in a sorted
#	list. The interval from the list is returned (or None).
#	Note that we consider abutting intervals as NON-overlapping.
#
# nota bene: this is a linear search; its performance could be improved by
#            using a binary search or a more sophisticated data structure

def find_overlap(intervals,start,end):
	# nota bene we rely on the fact that the intervals have been sorted
	for (s,e) in intervals:
		if (s >= end): return None
		if (e > start): return (s,e)
	return None


# count_overlap_bases--
#	Count the number of bases in an interval that overlap any interval in a
#	sorted list.
#
# nota bene: we assume the list of intervals does not contain any overlaps.

def count_overlap_bases(intervals,start,end):
	overlap = 0
	for (s,e) in intervals:
		if (e <= start): continue
		if (s >= end): break
		sOverlap = max(s,start)
		eOverlap = min(e,end)
		overlap += eOverlap - sOverlap
	return overlap


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


# parse_probability--
#	Parse a string as a probability

def parse_probability(s,strict=True):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	p *= scale

	if (strict) and (not 0.0 <= p <= 1.0):
		raise ValueError

	return p


if __name__ == "__main__": main()
