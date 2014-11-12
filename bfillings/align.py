# ----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from cogent import DNA as DNA_cogent, LoadSeqs
from cogent.align.align import make_dna_scoring_dict, local_pairwise

def pair_hmm_align_unaligned_seqs(seqs, moltype=DNA_cogent, params={}):
    """
        Checks parameters for pairwise alignment, returns alignment.

        Code from Greg Caporaso.
    """

    seqs = LoadSeqs(data=seqs, moltype=moltype, aligned=False)
    try:
        s1, s2 = seqs.values()
    except ValueError:
        raise ValueError(
            "Pairwise aligning of seqs requires exactly two seqs.")

    try:
        gap_open = params['gap_open']
    except KeyError:
        gap_open = 5
    try:
        gap_extend = params['gap_extend']
    except KeyError:
        gap_extend = 2
    try:
        score_matrix = params['score_matrix']
    except KeyError:
        score_matrix = make_dna_scoring_dict(
            match=1, transition=-1, transversion=-1)

    return local_pairwise(s1, s2, score_matrix, gap_open, gap_extend)
