import numpy as np
import os
import scrappy
import tempfile

from error_simulation import modelling

from fast5_research import Fast5
from uuid import uuid4


def simulate_errors(seq, basic=True):
    """
    Simulates errors occuring on given DNA sequence during the process of synthesis & sequencing.
    Parameters
    ----------
    seq : string
        DNA sequence errors will be simulated on
    basic : bool
        Determines whether to use simple simulation using predefined error rates (no quality scores),
        or a more advanced model which simulates signal generation and basecalling.
    Returns
    -------
    string
        DNA sequence containing simulated errors from synthesis and sequencing process.
    """
    print("Simulating synthesis errors...")
    syn_data = simulate_synthesis(seq)

    if basic:
        print("Simulating sequencing errors...")
        seq_data = simulate_sequencing(syn_data)
        return ''.join(seq_data)
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            print("Simulating sequencing...")
            simulate_read(syn_data, tmpdir, 'read.fast5')

            print("Simulating basecalling...")
            simulate_basecalling(tmpdir, tmpdir, 'basecalled.fastq')
            
            with open(os.path.join(tmpdir, 'basecalled.fastq')) as basecalled:
                lines = basecalled.readlines()
                out_seq, qscores = lines[1].strip(), lines[3].strip()
                return out_seq


def simulate_synthesis(seq, sub_rate=0.001, ins_rate=0.0015, del_rate=0.0055):
    """
    Simulates synthesis errors on sequence, outputs sequence with errors typical of synthesis
    process. Built upon existing modelling code.
    Parameters
    ----------
    seq : string
        List of nucleotide sequences to simulate synthesis errors on
    sub_rate : float, optional
        Substitution error rate for synthesis
    ins_rate : float, optional
        Insertion error rate for synthesis
    del_rate : float, optional
        Deletion error rate for synthesis
    Returns
    -------
    string
        List of nucleotide sequences with simulated synthesis errors
    """
    error_config = {
      'sub_rate': sub_rate,
      'ins_rate': ins_rate,
      'del_rate': del_rate
    }

    raw_data = modelling.decode_nucs(seq)
    syn_oligos, _ = modelling.synthesis(raw_data, error_config)
    return ''.join(modelling.encode_nucs(syn_oligos))


def simulate_sequencing(seq, sub_rate=0.15, ins_rate=0.05, del_rate=0.05):
    """
    Simulates sequencing errors on sequence, outputs sequence with errors typical of sequencing
    process. Built upon existing modelling code.
    Parameters
    ----------
    seq : string
        List of nucleotide sequences to simulate sequencing errors on
    sub_rate : float, optional
        Substitution error rate for sequencing
    ins_rate : float, optional
        Insertion error rate for sequencing
    del_rate : float, optional
        Deletion error rate for sequencing
    Returns
    -------
    string
        List of nucleotide sequences with simulated synthesis errors
    """
    error_config = {
      'sub_rate': sub_rate,
      'ins_rate': ins_rate,
      'del_rate': del_rate
    }

    raw_data = modelling.decode_nucs(seq)
    seq_oligos = modelling.sequence(raw_data, error_config)
    return ''.join(modelling.encode_nucs(seq_oligos))


def simulate_read(seq, out_dir, filename, id='test_read'):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    seq : string
        Nucleotide sequence
    out_dir : string
        Working directory to write intermediate results to.
    filename : string
        Name of FAST5 file which will be written as result of read.
        directory by default.
    id : string, optional
        ID given to read
    """

    # Using https://nanoporetech.github.io/fast5_research/examples.html as a reference
    squiggle = scrappy.sequence_to_squiggle(seq, rescale=True).data(as_numpy=True)
    raw_data = np.array([])

    for dwell, mean, stdv in squiggle:
        raw_data = np.append(raw_data, np.random.laplace(mean, stdv/np.sqrt(2), int(round(dwell))))

    start, stop = int(min(raw_data - 1)), int(max(raw_data + 1))
    rng = stop - start
    digitisation = 8192.0
    bins = np.arange(start, stop, rng / digitisation)
    # np.int16 is required, the library will refuse to write anything other
    raw_data = np.digitize(raw_data, bins).astype(np.int16)

    # The following are required meta data
    channel_id = {
        'digitisation': digitisation,
        'offset': 0,
        'range': rng,
        'sampling_rate': 4000,
        'channel_number': 1,
        }
    read_id = {
        'start_time': 0,
        'duration': len(raw_data),
        'read_number': 1,
        'start_mux': 1,
        'read_id': id,
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': id,
        'flow_cell_id': 'FAH00000',
    }
    context_tags = {}

    with Fast5.New(os.path.join(out_dir, filename), 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
        h.set_raw(raw_data, meta=read_id, read_number=1)


def simulate_basecalling(reads_dir, out_dir, out_file='basecalled.fastq'):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    reads_dir : string
        Directory from which all FAST5 files will be basecalled.
    out_file : string, optional
        File to write basecalled results to.
    """
    out_path = os.path.join(out_dir, out_file)
    os.system("bonito basecaller --fastq dna_r9.4.1@v2 %s > %s" % (reads_dir, out_path))
