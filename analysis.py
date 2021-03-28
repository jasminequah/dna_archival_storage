import pysam
import os

def count_errors(ref, alignment):
    """
    Computes alignment between reference sequence and sequence provided. Assumes performing single
    sequence alignment.
    Parameters
    ----------
    ref : string
        Name of fasta file containing reference sequence
    alignment : string
        Name of sam file containing alignment
    Returns
    -------
    (int, int, int, int):
        Number of insertions, deletions, substitutions, and length of reference.
    """
    with open(os.path.join(ref)) as ref:
        reference = ref.readlines()[1].strip()

    samfile = pysam.AlignmentFile(alignment, 'r')
    # Following https://samtools.github.io/hts-specs/SAMtags.pdf for parsing
    # CIGAR string. bwa-mem2 produces non-extended CIGAR string, so
    # substitution errors not seen in CIGAR string.
    for read in samfile.fetch(until_eof=True):
        alignment_start = read.reference_start
        alignment_end = read.reference_end
        ins = 0
        dels = alignment_start
        subs = 0
        read_pos = 0
        ref_pos = alignment_start

        for cigar_op, length in read.cigartuples:
            if cigar_op == 1:   # BAM_CINS (insertion)
                ins += length
                read_pos += length
            elif cigar_op == 2: # BAM_CDEL (deletion)
                dels += length
                ref_pos += length
            elif cigar_op == 4: # BAM_CSOFT_CLIP (soft clipping)
                if read_pos == 0:
                    dels -= length
                    ref_pos -= length
                subs += compute_subs(read_pos, ref_pos, length, reference, read.query_sequence)
                read_pos += length
                ref_pos += length
            else: # Assume BAM_CMATCH (alignment match)
                subs += compute_subs(read_pos, ref_pos, length, reference, read.query_sequence)
                read_pos += length
                ref_pos += length

        dels += len(reference) - ref_pos

        # md_tag = read.get_tag("MD")
        # subs += sum([1 for c in md_tag if c.isalpha()]) - dels
        # dels += alignment_start - 1 + reference_length - alignment_end

        print("ins: %d, dels: %d, subs: %d" % (ins, dels, subs))

    samfile.close()
    return ins, dels, subs, len(reference)


def compute_subs(read_pos, ref_pos, length, ref, read):
    subs = sum([1 for i in range(length) if ref[ref_pos + i] != read[read_pos + i]])
    return subs
