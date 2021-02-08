# dna_archival_storage

## Setup
Only standard libraries used (for now).

## `simulate_errors.py`
This script uses the `error_simulation` module to simulate synthesis and sequencing errors on a nucleotide sequence. It builds upon the work from a previous MEng student's project to model these errors. The libraries use the error profiles defined in their project report, which we will ideally replace with ours.

## Future Work (in progress)
* The DNA sequence on which errors are simulated is currently hardcoded, I will update the script to take a file as an argument
* Current error simulation assumes indels & substitution errors iid among sequence. Currently experimenting with trying to simulate errors using Scrappie for signal generation & Bonito basecalling.
