from encode import *
import filecmp


def test_bin_ascii_conversions():
    ascii_to_binary()
    binary_to_ascii()
    assert(filecmp.cmp('encoding_data/data.txt', 'encoding_data/output.txt'))


def test_nuc_conversions():
    ascii_to_binary()
    binary_to_nt()
    nt_to_binary()
    binary_to_ascii()
    assert(filecmp.cmp('encoding_data/data.txt', 'encoding_data/output.txt'))

def test_nuc_conversions_max_len_600():
    ascii_to_binary(max_len=600)
    binary_to_nt()
    nt_to_binary()
    binary_to_ascii()
    assert(filecmp.cmp('encoding_data/data.txt', 'encoding_data/output.txt'))
