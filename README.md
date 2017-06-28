# setMultiplex.py
Generates multiplexing designs for amplicon-sequencing diversity analysis using tagged primers

## Usage
```
python setMultiplex.py [-h] -i [I] -o [O] -f [F [F ...]] -r [R [R ...]] -n [N] [-l L L] [-m [M]] [-s [S]] [-rm [RM [RM ...]]] [-rmCombi [RMCOMBI [RMCOMBI ...]]] [-set [SET [SET ...]]] [--one] [--nolog] [--random] [--brute] [--verbose]
```

### Optional arguments
```
    -h, --help            Show help message and exit
    -i [FILE]             Tagged primers fasta file
    -o [FILE]             Output file name
    -f [F [F ...]]        Forward primer generic name
    -r [R [R ...]]        Reverse generic name
    -n [N]                Number of combinations to design. Enter different
                          numbers in the same order of the different tags files
    -l L L                Lengths of the forward and reverse tags
                          (default = [8, 8])
    -m [M]                Set a maximum primer usage frequency (automatically
                          changed if too low to be possible) (default = 0)
    -s [S]                Step-spacing between subsequent primers (if case of
                          sequencial cross-contamination of the primers and if
                          alphanumeric primer names) (default = 3)
    -rm [RM [RM ...]]     Single tagged primers not to include (either provide
                          file(s) with single primer name(s) in line(s), or
                          write it fully as space-separated arguents e.g. "V4F-4
                          F1-A")
    -rmCombi [RMCOMBI [RMCOMBI ...]]
                          File(s) with the combinations to avoid (fully written,
                          space separated)
    -set [SET [SET ...]]  Restict the choice of multiplexing design(s) to a list
                          of combinations (file(s) with forward and reverse
                          primers names in 1st and 2nd columns, respectively
    --one                 Get the first selection only (default: off)
    --nolog               Write a log file with the tested combination sets
                          (default: on)
    --random              Make random selection of mulitplexing design (default:
                          off)
    --brute               Make selection using brute-force algorithm based on
                          forward primers selection first and then reverse
                          primers until criteria satisfied (default: off)
    --verbose             Print ongoing processes and results to stdout
                          (default: off)
```

#### Requirements
Python2.7<br />
Biopython<br />
