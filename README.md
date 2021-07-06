# TrinoVap
->hier DOI

is a simple command line based tool to check for SwissProt, Kegg or Pfam annotations for a transcript based on the [Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki "Trinotate @ GitHub") annotation report.
TrinoVap also uses any [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki "TransDecoder @ GitHub") output in FASTA format to report annotions for complete, internal and partial sequences seperately.


## Requirements
[Python](https://www.python.org/downloads "Download Python") 2.7 or 3


## Usage
```bash
./TrinoVap.py TransDecoder FASTA Trinotate annotation report
```
### Positional Arguments:
```
  INPUT       name of TransDecoder FASTA (must be first argument),
              name of Trinotate annotation report TXT file (must be second argument)
```
