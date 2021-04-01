# dna_archival_storage

## Setup
To clone repository and submodules:
```
git clone --recursive https://github.com/jasminequah/dna_archival_storage
```

If running basic simulation (no basecalling), only standard libraries used.

If running simulation with raw signal simulation using Scrappie & Bonito basecalling, run `./setup.sh` first. A GPU is required to run the Bonito basecaller.


## Usage
Two modes of error simulation are currently available: one using given error rates and the `error_simulation` package (see below for more detail), and the other using Scrappie for signal generation and Bonito basecalling.

In Python, to run the simulation:
```
from simulation import simulate_errors
simulate_errors('AGGAATCTAGGCAGTAATAAATACATCAATCAATCAACTTAGCTATGCATTCATGAATAG', True)
```
Pass `True` for the second parameter to run the basic simulation (default), and `False` to run the simulation using Scrappie and Bonito.

## `error_simulation`
This package is used to simulate simple synthesis and sequencing errors on a nucleotide sequence for the basic simulation. It builds upon the work from a previous MEng student's project to model these errors. It uses the error profiles defined in their project report, which we will ideally replace with ours. It also assumes indels & substitution errors are iid among the sequence.

## Scrappie & Bonito simulation:
Simulates synthesis errors using same method as before, but simulates generated signals from sequencing process using Scrappie for R9.4.1 pore & basecalls these using Bonito basecaller with CTC architecture. Could also compare if DeepSimulator gives more representative simulation.

## TODOs
* Validate both types of error simulation with real data we have synthesised & sequenced
* Update error rates for basic simulation based on real data
* Maybe get the posterior probabilities for each possible base?
