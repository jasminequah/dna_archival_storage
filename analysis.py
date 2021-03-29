import pysam
import os

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
    (int, int, int, int):
        Number of insertions, deletions, substitutions, and length of reference.
    """
    if is_seq:
        # Create FASTA files for read and reference
        tmpdir = tempfile.TemporaryDirectory()
        with open(os.path.join(tmpdir, 'ref.fasta'), 'w') as ref_fasta:
            ref_fasta.write('>%s\n%s\n' % ('seq', ref))
        with open(os.path.join(tmpdir, 'read.fasta'), 'w') as read_fasta:
            read_fasta.write('>%s\n%s\n' % ('seq', read))
        ref = os.path.join(ref, 'ref.fasta')
        read = os.path.join(ref, 'read.fasta')

    alignment = align(ref, read)
    error_summary = count_errors(ref, alignment)

    if is_seq:
        tmpdir.cleanup()

    ins, dels, subs, ref_length = error_summary
    print("Insertion errors %.3f%%, deletion errors %.3f%%, substitution errors %.3f%%, error rate %.3f%%"
        % (ins / ref_length * 100, dels / ref_length * 100, subs / ref_length * 100, (ins + dels + subs) / ref_length * 100))

    # TODO: BLAST identity
    return error_summary


def align(ref, read, alignment="alignment.sam"):
    """
    Performs sequence alignment between ref and read using bwa-mem2.
    Parameters
    ----------
    ref : string
        Reference sequence filename
    read : string
        Read sequence filename
    alignment : string
        Name of SAM file which will be output, containing the sequence alignment.
    Returns
    -------
    string:
        Name of SAM file output.
    """
    os.system("./bwa-mem2-2.1_x64-linux/bwa-mem2 index %s" % ref)
    os.system("./bwa-mem2-2.1_x64-linux/bwa-mem2 mem %s %s > %s" % (ref, read, alignment))
    # TODO: Add setup instructions for bwa-mem2 to readme
    return alignment


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
        ins = 0
        dels = alignment_start
        subs = 0
        read_pos = 0
        ref_pos = alignment_start

        if read.cigartuples is None:
            raise Exception("Could not produce alignment")

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

        # Insertions and deletions at end of sequence
        dels += max(0, len(reference) - ref_pos)
        ins  += max(0, ref_pos - len(reference))

        # md_tag = read.get_tag("MD")
        # subs += sum([1 for c in md_tag if c.isalpha()]) - dels
        # dels += alignment_start - 1 + reference_length - alignment_end

        print("ins: %d, dels: %d, subs: %d" % (ins, dels, subs))

    samfile.close()
    return ins, dels, subs, len(reference)


def compute_subs(read_pos, ref_pos, length, ref, read):
    subs = sum([1 for i in range(length) if ref_pos + i < len(ref) and ref[ref_pos + i] != read[read_pos + i]])
    return subs
