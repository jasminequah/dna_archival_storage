import os
import numpy as np
import matplotlib.pyplot as plt

def get_error_summary(ref, read, is_seq=False):
    """
    Computes alignment between reference sequence and sequence provided. Assumes performing single
    sequence alignment.
    Parameters
    ----------
    ref : string
        Reference sequence
    read : string
        Read sequence
    is_seq : bool
        If True, ref and read are interpreted as raw sequences. Otherwise, they are interpreted as
        filenames for the reference and read sequence.
    Returns
    -------
    (ndarray(dtype=int), ndarray(dtype=int), ndarray(dtype=int), string, string):
        Positions of insertions, deletions, substitutions in sequence and alignment of ref and read.
    """
    if not is_seq:
        # Read from files
        with open(os.path.join(ref)) as ref_f:
            ref = ref_f.readlines()[1].strip()
        with open(os.path.join(read)) as read_f:
            read = read_f.readlines()[1].strip()

    error_summary = sm_align(ref, read)

    ref_align, read_align, _, ins_pos, dels_pos, subs_pos = error_summary
    ref_length = len(ref)
    ins_count = np.sum(ins_pos)
    dels_count = np.sum(dels_pos)
    subs_count = np.sum(subs_pos)
    print("Insertion errors %.3f%%, deletion errors %.3f%%, substitution errors %.3f%%, error rate %.3f%%"
        % (ins_count / ref_length * 100, dels_count / ref_length * 100, subs_count / ref_length * 100, (ins_count + dels_count + subs_count) / ref_length * 100))

    return ins_pos, dels_pos, subs_pos, ref_align, read_align


def get_errors(ref, read, is_seq=False):
    if len(read) == 0:
        print("LOST SEQUENCE: ignoring errors")
        return 0, 0, 0, 0
    ins_pos, dels_pos, subs_pos, ref_align, read_align = get_error_summary(ref, read, is_seq)
    ref_length = len(ref)
    ins_count = np.sum(ins_pos)
    dels_count = np.sum(dels_pos)
    subs_count = np.sum(subs_pos)
    return ins_count / ref_length * 100, dels_count / ref_length * 100, subs_count / ref_length * 100, (ins_count + dels_count + subs_count) / ref_length * 100


# get_similarity_score and sm_align adapted from
# https://github.com/JiaShun-Xiao/BLAST-bioinfor-tool/blob/master/blast.py
# TODO: refactor
def get_similarity_score(a, b, match_score=1, mismatch_score=-1):
    if a == b:
        return 1
    else:
        return -1


# Smith–Waterman Alignment
def sm_align(seq1, seq2, gap_penalty=3):
    m = len(seq1)
    n = len(seq2)
    ins_pos = np.zeros(m)
    dels_pos = np.zeros(m)
    subs_pos = np.zeros(m)
    g = -gap_penalty
    matrix = []
    for i in range(0, m):
        tmp = []
        for j in range(0, n):
            tmp.append(0)
        matrix.append(tmp)
    for sii in range(0, m):
        matrix[sii][0] = sii*g
    for sjj in range(0, n):
        matrix[0][sjj] = sjj*g
    for siii in range(1, m):
        for sjjj in range(1, n):
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + get_similarity_score(seq1[siii], seq2[sjjj]), matrix[siii][sjjj-1] + g)
    sequ1 = [seq1[m-1]]
    sequ2 = [seq2[n-1]]

    while m > 1 and n > 1:
        if max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-2][n-2]:
            m -= 1
            n -= 1
            sequ1.append(seq1[m-1])
            sequ2.append(seq2[n-1])
            if seq1[m-1] != seq2[n-1]:
                # Subsitution
                subs_pos[m-1] +=1
        elif max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-1][n-2]:
            # Insertion
            n -= 1
            sequ1.append('-')
            sequ2.append(seq2[n-1])
            ins_pos[m-2] += 1
        else:
            # Deletion
            m -= 1
            sequ1.append(seq1[m-1])
            sequ2.append('-')
            dels_pos[m-1] += 1
    sequ1.reverse()
    sequ2.reverse()
    align_seq1 = ''.join(sequ1)
    align_seq2 = ''.join(sequ2)
    align_score = 0.
    for k in range(0, len(align_seq1)):
        if align_seq1[k] == align_seq2[k]:
            align_score += 1
    align_score = float(align_score)/len(align_seq1)
    return align_seq1, align_seq2, align_score, ins_pos, dels_pos, subs_pos


def get_errors_from_file(filename, outdir='', limit=200):
    ins_pos = np.array([]) # relative
    dels_pos = np.array([]) # relative
    subs_pos = np.array([]) # relative
    ins_rates = []
    dels_rates = []
    subs_rates = []
    error_rates = []
    ref_lengths = []
    read_lengths = [] # relative
    with open(filename, 'r') as f:
        num_samples = 0
        line = f.readline()
        while line and num_samples < limit:
            ref = f.readline().strip()
            read = f.readline().strip()
            error_summary = sm_align(ref, read)
            ref_align, read_align, _, ins, dels, subs = error_summary
            with open(os.path.join(outdir, 'error_pos.log'), 'a') as log:
                log.write("%s\n%s\n%s\n" % (ins.tolist(), dels.tolist(), subs.tolist()))

            relative_pos = lambda x: x / len(ref)
            ins_rates.append(np.sum(ins) / len(ref) * 100)
            ins = [[i for j in range(int(err_count))] for i, err_count in enumerate(ins)]
            ins = [relative_pos(pos) for sublist in ins for pos in sublist]
            ins_pos = np.concatenate((ins_pos, ins))

            dels_rates.append(np.sum(dels) / len(ref) * 100)
            dels = [[i for j in range(int(err_count))] for i, err_count in enumerate(dels)]
            dels = [relative_pos(pos) for sublist in dels for pos in sublist]
            dels_pos = np.concatenate((dels_pos, dels))

            subs_rates.append(np.sum(subs) / len(ref) * 100)
            subs = [[i for j in range(int(err_count))] for i, err_count in enumerate(subs)]
            subs = [relative_pos(pos) for sublist in subs for pos in sublist]
            subs_pos = np.concatenate((subs_pos, subs))

            error_rates.append(ins_rates[-1] + dels_rates[-1] * subs_rates[-1])
            print("%s, %s, %s, %s" % (ins_rates[-1], dels_rates[-1], subs_rates[-1], ins_rates[-1] + dels_rates[-1] + subs_rates[-1]))

            ref_lengths.append(len(ref))
            read_lengths.append(len(read) / len(ref))

            line = f.readline()
            num_samples += 1

    ins_rate = np.mean(ins_rates)
    dels_rate = np.mean(dels_rates)
    subs_rate = np.mean(subs_rates)
    print("Average insertions %.3f%%, deletions %.3f%%, substitutions %.3f%%, errors %.3f%%"
        % (ins_rate, dels_rate, subs_rate, ins_rate + dels_rate + subs_rate))

    # Plotting distribution of insertion, deletion & substitution errors
    plt.figure()
    plt.title("Sequence Error Distribution")
    x, y = np.unique(ins_pos, return_counts=True)
    x, y = bucket_frequencies(x, y, 0, 1, 0.002)
    plt.plot(x, y, label="Insertions")
    x, y = np.unique(dels_pos, return_counts=True)
    x, y = bucket_frequencies(x, y, 0, 1, 0.002)
    plt.plot(x, y, label="Deletions")
    x, y = np.unique(subs_pos, return_counts=True)
    x, y = bucket_frequencies(x, y, 0, 1, 0.002)
    plt.plot(x, y, label="Substitutions")
    plt.xlabel('DNA Sequence Relative Position')
    plt.ylabel('Total Error Count')
    plt.legend()
    plt.savefig(os.path.join(outdir, "error_distribution.png"))

    # Write statistics
    with open(os.path.join(outdir, 'avg_error_summary.txt'), 'w') as fasta:
        fasta.write("Insertions(%%):\t\t%.3f\nDeletions(%%):\t\t\t%.3f\nSubstitutions(%%):\t%.3f\nErrors(%%):\t\t\t\t%.3f\n"
            % (ins_rate, dels_rate, subs_rate, (ins_rate + dels_rate + subs_rate)))

    # Plot read length vs error rate
    plt.figure()
    plt.title("Sequence Length vs Error Rate")
    plt.plot(ref_lengths, error_rates, 'o')
    plt.plot(np.unique(ref_lengths), np.poly1d(np.polyfit(ref_lengths, error_rates, 1))(np.unique(ref_lengths)))
    plt.xlabel('Sequence Length (nt)')
    plt.ylim(0,10)
    plt.ylabel('Error Rate (%)')
    plt.savefig(os.path.join(outdir, "len_vs_errors.png"))

    # Plot distribution of read lengths
    plt.figure()
    x, y = np.unique(read_lengths, return_counts=True)
    x, y = bucket_frequencies(x, y, 0.88, 1.12, 0.005)
    x, y = filter_zero_counts(x, y)
    plt.title("Relative Read Length Distribution")
    plt.scatter(x, y)
    plt.ylabel('Frequency')
    plt.xlabel('Relative Read Length')
    plt.axvline(x=1, label="Reference Length")
    plt.savefig(os.path.join(outdir, "read_len_distribution.png"))


def bucket_frequencies(x, y, lo, hi, interval):
    x_new = np.arange(lo + interval / 2, hi, interval)
    y_new = np.zeros(int((hi - lo) / interval))
    for i in range(len(x)):
        y_new[int((x[i] - lo) // interval)] += y[i]
    return x_new, y_new


def filter_zero_counts(x, y):
    x_new = []
    y_new = []
    for i in range(len(y)):
        if y[i] != 0:
            x_new.append(x[i])
            y_new.append(y[i])
    return x_new, y_new
