# fineRADstructure-tools
Tools for data conversion and results visualization for fineRADstructure for fineRADstructure (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html)

## Tool 1: finerad_input.py
Type `python finerad_input.py -h` to show the help:
```
usage: finerad_input.py [-h] -i FILENAME [-t DATA_TYPE] [-n MINSAMPLE]
                        [-o OUTFILE]

Converts haplotype data from Stacks, pyrad or ipyrad into a haplotype matrix
for analysis with fineRADstructure
(http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html)

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the haplotypes.tsv file from Stacks, the
                        .alleles file from pyrad, or the .alleles.loci file
                        from ipyrad to be converted.
  -t DATA_TYPE, --type DATA_TYPE
                        Source of the input file, options: stacks, pyrad,
                        ipyrad, default=guess
  -n MINSAMPLE, --minsample MINSAMPLE
                        Minimum number of samples in a locus, default=2
  -o OUTFILE, --output OUTFILE
                        Name for the output file, default=input name +
                        .minsample + .finerad
```

_Example 1:_ Convert the .alleles.loci matrix from ipyrad and exclude loci with less than 8 samples:
```
python finerad_input.py --input example_data.alleles.loci --minsample 8
984 loci for 12 samples written to example_data.alleles.loci.min8.finerad
```

_Example 2:_ Convert the .haplotypes.tsv matrix from Stacks, exclude loci with 4 samples or less and specify name of output matrix:
```
python finerad_input.py -i batch_1.haplotypes.tsv -n 4 --output stacks_finerad_min4.finerad
152181 loci for 37 samples written to stacks_finerad_min4.finerad
```

_Example 3:_ Convert the .alleles matrix from pyrad with default options:
```
python finerad_input.py -i pyrad_data.alleles
13007 loci for 170 samples written to pyrad_data.alleles.min2.finerad
```

## References
**Malinsky, M., E. Trucchi, D. Lawson & D. Falush. 2016.** RADpainter and fineRADstructure: population inference from RADseq data. _bioRxiv_. doi:https://doi.org/10.1101/057711

A description of `fineRADstructure` and links to other relevant papers is available at: http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html  
`fineRADstructure` code is available at: https://github.com/millanek/fineRADstructure

## Credits
- Code: [Edgardo M. Ortiz](mailto:e.ortiz.v@gmail.com)

