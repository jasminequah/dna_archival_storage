# import scrappy
# import numpy as np
# import bonito

from error_simulation import modelling

# from uuid import uuid4
# from fast5_research import Fast5

# Encoder input: array of letters in {A, C, T, G} & list determining mapping from {0,1}^2 to A,C,T,G
# Encoder output: array of letters from {A, C, T, G} which may be longer

# Decoder input: array of letters in {A, C, T, G, E} & confidence
# Decoder output: array of letters in {A, C, T, G}

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
    string[]
        List of nucleotide sequences with simulated synthesis errors
    """
    error_config = {
      'sub_rate': sub_rate,
      'ins_rate': ins_rate,
      'del_rate': del_rate
    }

    raw_data = modelling.decode_nucs(seq)
    syn_oligos, _ = modelling.synthesis(raw_data, error_config)
    return modelling.encode_nucs(syn_oligos)


def simulate_sequencing(seq, sub_rate=0.15, ins_rate=0.05, del_rate=0.05):
    """
    Simulates synthesis errors on sequence, outputs sequence with errors typical of sequencing
    process. Built upon existing modelling code.
    Parameters
    ----------
    seq : string[]
        List of nucleotide sequences to simulate sequencing errors on
    sub_rate : float, optional
        Substitution error rate for sequencing
    ins_rate : float, optional
        Insertion error rate for sequencing
    del_rate : float, optional
        Deletion error rate for sequencing
    Returns
    -------
    string[]
        List of nucleotide sequences with simulated synthesis errors
    """
    error_config = {
      'sub_rate': sub_rate,
      'ins_rate': ins_rate,
      'del_rate': del_rate
    }

    raw_data = modelling.decode_nucs(seq)
    seq_oligos = modelling.sequence(raw_data, error_config)
    return modelling.encode_nucs(seq_oligos)


def simulate_read(seq, filename='reads/test.fast5'):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    seq : string
        Nucleotide sequence
    filename: string
        Name of FAST5 file which will be written as result of read
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
        'read_id': 'test_read',
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': 'test_read',
        'flow_cell_id': 'FAH00000',
    }
    context_tags = {}

    with Fast5.New(filename, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
        h.set_raw(raw_data, meta=read_id, read_number=1)


def simulate_basecalling(file):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    file : string
        Name of FAST5 file input to basecaller
    Returns
    -------
    TODO
    """
    pass


def main():

  data = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
  print("Encoded data:\n %s" % data)

  print("Simulating synthesis errors...")
  syn_data = simulate_synthesis(data)

  print("Synthesised data:\n %s" % ''.join(syn_data))

  print("Simulating sequencing errors...")
  seq_data = simulate_sequencing(syn_data)

  print("Data read back:\n %s" % ''.join(seq_data))


if __name__ == "__main__":
  main()
