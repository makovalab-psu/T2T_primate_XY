# Palindrover

## System Requirements

This requires Python 3 and has been tested on a common Linux environment.

## Installation

If you have Python 3, Palindrover should work. Installation should be as fast
as downloading this repository. Start by downloading the repository:

```
git clone https://github.com/makovalab-psu/T2T_primate_XY
cd T2T_primate_xy/src
```

Simply run by providing the path to `palindrover.py`:

```
cat /path/to/self-alignments.lz | python3 ./palindrover.py [options]
```

## Theoretical Example

Generate self-alignments with output in LASTZ format with cigarx. Alignemnts
can be from any aligner as long as th format is correct. An example from LASTZ
is as follows:

```
lastz input.fa --self \
	--format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
	> self-alignments.lz
```

Then, run palindrover (runs in seconds):

```
python3 ./palindrover.py  \
	< self-alignments.lz \
	> self-alignments.palindromes.lz
```

## How to run for ape XY data

Running it on the ape XY data is the same as the example, simply provide an
alignment file in LASTZ format on standard in. It would need to be done on
self-alignments for each species for each chromosome. Here's an example for
Chimpanzee chrY:

```
lastz chimp.chrY.fa --self \
	--format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
	> chimp-chrY-self-alignments.lz

python3 ./palindrover.py  \
	< chimp-chrY-self-alignments.lz \
	> chimp-chrY-self-alignments.palindromes.lz
```

## Testing

You can test your installation on a test dataset. A subset of the bonobo chrY
self-alignments from LASTZ are included in the `test` directory. Try running
like this (runs in seconds):

```
python3 ./palidnrover.py \
	< test/bonoboY-selfAlns.lz \
	> bonoboY-selfAlns.palindromes.lz
```

Then check if your output matches the expected output:

```
diff -sq bonoboY-selfAlns.palindromes.lz test/bonoboY-selfAlns.palindromes.lz
```

You should see the following output:

```
Files bonoboY-selfAlns.palindromes.lz and test/bonoboY-selfAlns.palindromes.lz are identical
```

Optionally, delete the test output you created:

```
rm bonoboY-selfAlns.palindromes.lz
```

