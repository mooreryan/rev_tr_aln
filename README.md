# Teeheehee #

## Overview ##

Given an AA alignment and a file with AA and NT (DNA) orfs, convert
the AA alignment file into a NT alignment.

## Dependencies ##

- `parse_fasta`
- `trollop`

### Install gems ###

    gem install parse_fasta
    gem install trollop

## Linking the ORFs ##

We assume that the headers in the AA alingment file match those in the
AA orfs file.

Headers for the AA orfs and NT orfs are broken on spaces and matched
on the first of these tokens. E.g.,

    >gi|456794300|gb|AMFO01000279.1|_18 blah blah

would match on `gi|456794300|gb|AMFO01000279.1|_18`, and

    >3242856 tee hee hee

would match on `3242856`.

## Start and stop codons ##

Some AA orfs and NT orfs don't match regarding the start codon. If one
has a stop and the other doesn't the next position is checked for
plausibility.

## Ambigious positions ##

- Amino acid `X` becomes `NNN`
- If the NT orf has an ambiguous base it should be skipped

## Codon table ##

Link to [codon table](http://www.cbs.dtu.dk/courses/27619/codon.html)
