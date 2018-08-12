## Augurlinos: a collection of modules for molecular epidemiology

_This repo was a prototype of how to strip down augur into a series of interacting modules. All development has moved to the [augur](https://github.com/nextstrain/augur) repository._

#### Quickstart
Put a `zika.fasta` into `zika/data` and run

```
snakemake --snakefile zika/zika_snake export
```

### Requirments
  * python 2.7
  * numpy, biopython, scipy
  * treetime

### Structure
Each build should have a Snakefile. The example zika/zika_snake requires a fasta file with zika sequence (fauna output) and an annotated reference sequence in genbank format. Each feature that should be translated needs a "locus_tag" which will be used to name the translation.

Different steps of a build produce standard files (alignments in fasta format, tree in newick format, tables as tsv) on which subsequent steps rely. The final step exports a jsons for auspice

### ToDos
  * meta.json is not yet created
  * frequencies and titers at not yet done (easy)
  * support for multiple segments is missing
  * prepare doesn't do anything meaningful.
