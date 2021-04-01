import numpy as np
import os
import random
import scrappy
import tempfile

from error_simulation import modelling
from analysis import get_error_summary

from fast5_research import Fast5
from uuid import uuid4


def simulate_errors(seq, basic=True, save_dir=None, id='test_read', print_error_summary=False):
    """
    Simulates errors occuring on given DNA sequence during the process of synthesis & sequencing.
    Parameters
    ----------
    seq : string
        DNA sequence errors will be simulated on
    basic : bool, optional
        Determines whether to use simple simulation using predefined error rates (no quality scores),
        or a more advanced model which simulates signal generation and basecalling.
    save_dir : string, optional
        File path to save intermediate fasta, fast5, and fastq files generated by simulation process
        into. If None, these are not saved
    id : string, optional
        ID given to read
    Returns
    -------
    string, [float], ?tuple
        DNA sequence containing simulated errors from synthesis and sequencing process,
        qscores associated with each, and error summary if print_error_summary is True.
    """
    print("Simulating synthesis errors...")
    syn_data = simulate_synthesis(seq)

    if basic:
        print("Simulating sequencing errors...")
        seq_data = simulate_sequencing(syn_data)
        return ''.join(seq_data), [0.25] * len(seq_data)
    else:
        working_dir = save_dir
        if working_dir is None:
            tmpdir = tempfile.TemporaryDirectory()
            working_dir = tmpdir.name

        with open(os.path.join(working_dir, 'seq.fasta'), 'w') as fasta:
            fasta.write('>%s\n%s\n' % (id, seq))

        print("Simulating sequencing...")
        simulate_read(syn_data, working_dir, 'read.fast5', id)

        print("Simulating basecalling...")
        simulate_basecalling(working_dir, working_dir, 'basecalled.fastq')

        with open(os.path.join(working_dir, 'basecalled.fastq')) as basecalled:
            lines = basecalled.readlines()
            out_seq, qscores = lines[1].strip(), lines[3].strip()

        # Print insertion, deletion, substitution error statistics
        error_summary = None
        if print_error_summary:
            print(seq)
            print(out_seq)
            error_summary = get_error_summary(os.path.join(working_dir, 'seq.fasta'), os.path.join(working_dir, 'basecalled.fastq'))

        if save_dir is None:
            tmpdir.cleanup()
        return out_seq, parse_qscores(qscores), error_summary


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


def simulate_read(seq, out_dir, filename, id):
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
    id : string
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


def parse_qscores(qscores):
    """
    Converts quality score character encoding into probability values
    """
    return [1 - 10 ** (-(ord(qscore) - 33)/10) for qscore in qscores]


def evaluate_simulator(seq_length=100, runs=20):
    avg_ins = 0
    avg_dels = 0
    avg_subs = 0
    for i in range(runs):
        seq = create_random_seq(seq_length)
        (_, _, error_summary) = simulate_errors(seq, False, print_error_summary=True)
        (ins, dels, subs, _) = error_summary
        avg_ins += ins
        avg_dels += dels
        avg_subs += subs
    avg_ins /= runs
    avg_dels /= runs
    avg_subs /= runs
    print("Average insertions %.3f%%, deletions %.3f%%, substitutions %.3f%%, errors %.3f%%"
        % (avg_ins / seq_length * 100, avg_dels / seq_length * 100, avg_subs / seq_length * 100, (avg_ins + avg_dels + avg_subs) / seq_length * 100))


def create_random_seq(seq_length):
    return ''.join(random.choices(['A','C','G','T'], k=seq_length))
