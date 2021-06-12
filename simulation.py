import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
import scrappy
import tempfile

from error_simulation import modelling
from analysis import get_error_summary

from fast5_research import Fast5
from uuid import uuid4


def simulate_errors(seqs, basic=True, save_dir=None, id='test_read', print_error_summary=False, error_rate=None, no_indels=False):
    """
    Simulates errors occuring on given DNA sequence during the process of synthesis & sequencing.
    Parameters
    ----------
    seqs : string or [string]
        DNA sequence(s) errors will be simulated on
    basic : bool, optional
        Determines whether to use simple simulation using predefined error rates (no quality scores),
        or a more advanced model which simulates signal generation and basecalling.
    save_dir : string, optional
        File path to save intermediate fasta, fast5, and fastq files generated by simulation process
        into. If None, these are not saved.
    id : string, optional
        ID prefix given to reads.
    print_error_summary : bool, optional
        If true and basic=False, prints error statistics of simulated sequence.
    error_rate : float, optional
        If set, either artifically increases or decreases the error rate towards the specified rate.
        Must be in range [0,1].
    no_indels : bool, optional
        If true, manually fixes all indel errors on simulated sequences, leaving only substitution errors.
    Returns
    -------
    [string], [[float]], [?tuple]
        DNA sequence(s) containing simulated errors from synthesis and sequencing process,
        qscores associated with each, and error summary if print_error_summary is True.
    """
    if not isinstance(seqs, list):
        seqs = [seqs]

    out_seqs = []
    out_qscores = []
    error_summaries = [None] * len(seqs)
    working_dir = save_dir
    tmpdir = None

    print("Simulating synthesis & sequencing...")
    for i, seq in enumerate(seqs):
        syn_data = simulate_synthesis(seq)

        if basic:
            seq_data, _ = simulate_sequencing(syn_data)
            out_seqs.append(''.join(seq_data))
            out_qscores.append([0.25] * len(seq_data))
        else:
            if working_dir is None:
                tmpdir = tempfile.TemporaryDirectory()
                working_dir = tmpdir.name

            with open(os.path.join(working_dir, 'seq_%s.fasta' % i), 'w') as fasta:
                fasta.write('>%s\n%s\n' % ("%s_%s" % (id, i), seq))

            simulate_read(syn_data, working_dir, 'read_%s.fast5' % i, "%s_%s" % (id, i))

    if basic:
        return out_seqs, out_qscores, error_summaries

    print("Simulating basecalling...")
    simulate_basecalling(working_dir, working_dir, 'basecalled.fastq')

    out_seqs = [None] * len(seqs)
    out_qscores = [None] * len(seqs)
    with open(os.path.join(working_dir, 'basecalled.fastq')) as basecalled:
        line = basecalled.readline()
        while line:
            read_index = int(re.match('^.*_([0-9]+)$', line.strip()).group(1))
            out_seqs[read_index] = basecalled.readline().strip()
            basecalled.readline()
            out_qscores[read_index] = basecalled.readline().strip()
            line = basecalled.readline()

    for i, seq in enumerate(seqs):
        if no_indels:
            out_seq, qscores = fix_indels(seq, out_seqs[i], out_qscores[i])
            out_seqs[i] = out_seq
            out_qscores[i] = qscores

        if error_rate is not None:
            out_seq, qscores = alter_error_rate(seq, out_seqs[i], out_qscores[i], error_rate, no_indels)
            out_seqs[i] = out_seq
            out_qscores[i] = qscores

        # Print insertion, deletion, substitution error statistics
        if print_error_summary:
            print(seq)
            print(out_seqs[i])
            error_summary = get_error_summary(seq, out_seqs[i], True)
            error_summaries[i] = error_summary

    if no_indels or error_rate is not None:
        # Update fastq file
        with open(os.path.join(working_dir, 'basecalled.fastq'), 'r') as fastq:
            lines = fastq.readlines()
        with open(os.path.join(working_dir, 'basecalled.fastq'), 'w') as fastq:
            for i, seq in enumerate(seqs):
                lines[i * 4 + 1] = out_seqs[i]
                lines[i * 4 + 3] = out_qscores[i]
            fastq.writelines(lines)

    if save_dir is None:
        tmpdir.cleanup()

    return out_seqs, [parse_qscores(qscores) for qscores in out_qscores], error_summaries


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
    string, (ndarray(dtype=bool), ndarray(dtype=bool), ndarray(dtype=bool))
        List of nucleotide sequences with simulated sequencing errors and tuple showing positions in
        original sequence where insertions, deletions, and substitutions occur.
    """
    error_config = {
        'sub_rate': sub_rate,
        'ins_rate': ins_rate,
        'del_rate': del_rate
    }

    raw_data = modelling.decode_nucs(seq)
    seq_oligos, errors_pos = modelling.sequence(raw_data, error_config)
    return ''.join(modelling.encode_nucs(seq_oligos)), errors_pos


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
    os.system("bonito basecaller --fastq --weights 1 dna_r9.4.1@v3.2 %s > %s" % (reads_dir, out_path))


def parse_qscores(qscores):
    """
    Converts quality score character encoding into probability values
    """
    return [1 - 10 ** (-(ord(qscore) - 33)/10) for qscore in qscores]


def evaluate_simulator(seq_length=100, runs=20, outdir='', error_rate=None, no_indels=False):
    ins_pos = np.zeros(seq_length)
    dels_pos = np.zeros(seq_length)
    subs_pos = np.zeros(seq_length)

    seqs = [create_random_seq(seq_length) for i in range(runs)]
    reads, _, error_summaries = simulate_errors(seqs, False, print_error_summary=True, error_rate=error_rate, no_indels=no_indels)
    read_lengths = [len(read) for read in reads]

    for i, error_summary in enumerate(error_summaries):
        (ins, dels, subs, ref_align, read_align) = error_summary
        np.add(ins_pos, ins, out=ins_pos)
        np.add(dels_pos, dels, out=dels_pos)
        np.add(subs_pos, subs, out=subs_pos)

        with open(os.path.join(outdir, 'simulated_seqs_debug.txt'), 'a+') as simulated_seqs:
            simulated_seqs.write(">simulated_seq_%s\n%s\n%s\n" % (i, seqs[i], reads[i]))
        with open(os.path.join(outdir, 'simulated_seqs_debug_align.txt'), 'a+') as simulated_seqs:
            ref_length = seq_length
            ins_count = np.sum(ins)
            dels_count = np.sum(dels)
            subs_count = np.sum(subs)
            error_txt = "Insertion errors %.3f%%, deletion errors %.3f%%, substitution errors %.3f%%, error rate %.3f%%" \
                % (ins_count / ref_length * 100, dels_count / ref_length * 100, subs_count / ref_length * 100, (ins_count + dels_count + subs_count) / ref_length * 100)
            simulated_seqs.write(">simulated_seq_%s\n%s\n%s\n%s\n" % (i, ref_align, read_align, error_txt))

    ins_rate = np.sum(ins_pos) / runs / seq_length * 100
    dels_rate = np.sum(dels_pos) / runs / seq_length * 100
    subs_rate = np.sum(subs_pos) / runs / seq_length * 100
    print("Average insertions %.3f%%, deletions %.3f%%, substitutions %.3f%%, errors %.3f%%"
        % (ins_rate, dels_rate, subs_rate, ins_rate + dels_rate + subs_rate))

    # Plotting distribution of insertion, deletion & substitution errors
    plt.figure()
    plt.title("Sequence Error Distribution")
    plt.plot(ins_pos, label="Insertions")
    plt.plot(dels_pos, label="Deletions")
    plt.plot(subs_pos, label="Substitutions")
    plt.xlabel('DNA Sequence Position')
    plt.ylabel('Total Error Count')
    plt.legend()
    plt.savefig(os.path.join(outdir, "simulated_error_distribution.png"))

    # Write statistics
    with open(os.path.join(outdir, 'simulated_avg_error_summary.txt'), 'w') as fasta:
        fasta.write("Insertions(%%):\t\t%.3f\nDeletions(%%):\t\t\t%.3f\nSubstitutions(%%):\t%.3f\nErrors(%%):\t\t\t\t%.3f\n"
            % (ins_rate, dels_rate, subs_rate, (ins_rate + dels_rate + subs_rate)))

    # Plot distribution of read lengths
    plt.figure()
    x, y = np.unique(read_lengths, return_counts=True)
    plt.title("Read Length Distribution")
    plt.scatter(x, y)
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.axvline(x=seq_length, label="Reference length")
    plt.savefig(os.path.join(outdir, "simulated_read_len_distribution.png"))


def create_random_seq(seq_length):
    return ''.join(random.choices(['A','C','G','T'], k=seq_length))


def fix_indels(ref, seq, qscores=None):
    print("Before fixing indels:")
    error_summary = get_error_summary(ref, seq, True)
    ins_pos, dels_pos, subs_pos, ref_align, read_align = error_summary
    align_len = len(ref_align)
    seq_pos = 0
    new_seq = []
    new_qscores = [] if qscores is not None else None
    for align_pos in range(align_len):
        if ref_align[align_pos] == '-':
            # Insertion
            seq_pos += 1
        elif read_align[align_pos] == '-':
            # Deletion
            new_seq.append(ref_align[align_pos])
            if new_qscores is not None:
                new_qscores.append('=') # TODO: pick qscore for correcting deletion
        else:
            new_seq.append(seq[seq_pos])
            if new_qscores is not None:
                new_qscores.append(qscores[seq_pos])
            seq_pos += 1
    new_qscores = ''.join(new_qscores) if qscores is not None else None
    return ''.join(new_seq), new_qscores


def alter_error_rate(ref, seq, qscores, desired_error_rate, no_indels=False):
    print("Before tuning error rate:")
    error_summary = get_error_summary(ref, seq, True)
    ins_pos, dels_pos, subs_pos, ref_align, read_align = error_summary
    curr_error_rate = (np.sum(ins_pos) + np.sum(dels_pos) + np.sum(subs_pos)) / len(ref)

    if curr_error_rate > desired_error_rate:
        # Correct a proportion of errors in sequence
        # Each error is corrected with probability (curr_error_rate - desired_error_rate) / curr_error_rate
        correction_prob = (curr_error_rate - desired_error_rate) / curr_error_rate
        align_len = len(ref_align)
        seq_pos = 0
        new_seq = []
        new_qscores = []
        for align_pos in range(align_len):
            # Correct error at alignment position
            if ref_align[align_pos] == '-':
                if random.random() < correction_prob:
                    # Insertion
                    if not no_indels:
                        seq_pos += 1
                        continue
            elif read_align[align_pos] == '-':
                if random.random() < correction_prob:
                    # Deletion
                    if not no_indels:
                        new_seq.append(ref_align[align_pos])
                        new_qscores.append('=') # TODO: pick qscore for correcting deletion
                        continue
            elif ref_align[align_pos] != read_align[align_pos]:
                if random.random() < correction_prob:
                    # Substitution
                    new_seq.append(ref_align[align_pos])
                    new_qscores.append(qscores[seq_pos]) # TODO: pick qscore for corrected base
                    seq_pos += 1
                    continue
            if read_align[align_pos] != '-':
                new_seq.append(seq[seq_pos])
                new_qscores.append(qscores[seq_pos])
                seq_pos += 1
        return ''.join(new_seq), ''.join(new_qscores)
    else:
        # Introduce errors with probability desired_error_rate - curr_error_rate. Assumes errors
        # uniformly distributed throughout sequence - could be improved in future versions.
        if no_indels:
            subs_rate = desired_error_rate - curr_error_rate
            ins_rate = 0.
            dels_rate = 0.
        else:
            subs_rate = (np.sum(subs_pos) / len(ref) / curr_error_rate) * (desired_error_rate - curr_error_rate)
            ins_rate = (np.sum(ins_pos) / len(ref) / curr_error_rate) * (desired_error_rate - curr_error_rate)
            dels_rate = (np.sum(dels_pos) / len(ref) / curr_error_rate) * (desired_error_rate - curr_error_rate)
        new_seq, errors_pos = simulate_sequencing(seq, subs_rate, ins_rate, dels_rate)

        # Updating qscores
        ins_pos_new, dels_pos_new, subs_pos_new = errors_pos
        new_qscores = []
        for pos in range(len(seq)):
            if ins_pos_new[pos]:
                new_qscores.append('=') # TODO: pick qscore for introducing insertion
            elif dels_pos_new[pos]:
                continue
            else:
                new_qscores.append(qscores[pos]) # TODO: assumes same qscore for substitution

        return new_seq, ''.join(new_qscores)
