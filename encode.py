def ascii_to_binary(f='encoding_data/data.txt', b='encoding_data/data_bin.txt', max_len=None):
    with open(f, 'r') as txtfile:
        text = txtfile.read()

    bin_seq = ''.join('{0:08b}'.format(ch, 'b') for ch in bytearray(text, 'utf-8'))
    if max_len:
        bin_seq = '\n'.join([bin_seq[0+i:max_len+i] for i in range(0, len(bin_seq), max_len)])

    with open(b, 'w') as binfile:
        binfile.write(bin_seq)

    return bin_seq


def binary_to_ascii(b='encoding_data/data_bin.txt', output='encoding_data/output.txt'):
    with open(b, 'r') as binfile:
        bin_seq = binfile.read()
    bin_seq = bin_seq.replace('\n', '')

    with open(output, 'w') as out:
        # Assuming utf-8
        decoded = ''.join(chr(int(bin_seq[i*8:i*8+8],2)) for i in range(len(bin_seq)//8))
        out.write(decoded)

    return decoded


def binary_to_nt(b='encoding_data/data_bin.txt', output='encoding_data/data_nuc.txt'):
    with open(b, 'r') as binfile:
        bin_seqs = binfile.readlines()

    dna_seqs = [binary_to_bases(bin_seq.strip()) for bin_seq in bin_seqs]

    with open(output, 'w') as out:
        out.write('\n'.join(dna_seqs))

    return dna_seqs


def nt_to_binary(nt='encoding_data/data_nuc.txt', output='encoding_data/data_bin.txt'):
    with open(nt, 'r') as ntfile:
        dna_seqs = ntfile.readlines()

    bin_seqs = [bases_to_binary(dna_seq.strip()) for dna_seq in dna_seqs]

    with open(output, 'w') as out:
        out.write('\n'.join(bin_seqs))

    return bin_seqs


def binary_to_bases(bin_seq):
    mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    return ''.join(mapping[bin_seq[i*2:i*2+2]] for i in range(len(bin_seq)//2))


def bases_to_binary(dna_seq):
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return ''.join(mapping[nt] for nt in dna_seq)
