# ----------------------------------------------------------------------------
# Copyright (c) 2015--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
"""
Provides an application controller for the commandline version of:
MAFFT v7.220
"""

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath
from skbio import Alignment, DNA, RNA, Protein

MOLTYPE_MAP = {DNA: '--nuc',
               RNA: '--nuc',
               Protein: '--amino'}


class Mafft(CommandLineApplication):
    """Mafft application controller"""
    _options = {
        # Algorithm

        # Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i
        # and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        '--auto': FlagParameter(Prefix='--', Name='auto'),

        # Distance is calculated based on the number of shared 6mers.
        # Default: on
        '--6merpair': FlagParameter(Prefix='--', Name='6merpair'),

        # All pairwise alignments are computed with Needleman-Wunsch algorithm.
        # More accurate but slower than --6merpair. Suitable for globally
        # alignable sequences. Applicable to up to ~200 sequences.
        # A combination with --maxiterate 1000 is recommended (G-INS-i).
        # Default: off (6mer distance is used)
        '--globalpair': FlagParameter(Prefix='--', Name='globalpair'),

        # All pairwise alignments are computed with the Smith-Waterman algorithm.
        # More accurate but slower than --6merpair. Suitable for a set of locally
        # alignable sequences. Applicable to up to ~200 sequences. A combination
        # with --maxiterate 1000 is recommended (L-INS-i). Default: off
        # (6mer distance is used)
        '--localpair': FlagParameter(Prefix='--', Name='localpair'),

        # All pairwise alignments are computed with a local algorithm with the
        # generalized affine gap cost (Altschul 1998). More accurate but slower than
        # --6merpair. Suitable when large internal gaps are expected. Applicable to
        # up to ~200 sequences. A combination with --maxiterate 1000 is recommended
        # (E-INS-i). Default: off (6mer distance is used)
        '--genafpair': FlagParameter(Prefix='--', Name='genafpair'),

        # All pairwise alignments are computed with FASTA (Pearson and Lipman 1988).
        # FASTA is required. Default: off (6mer distance is used)
        '--fastapair': FlagParameter(Prefix='--', Name='fastapair'),

        # Weighting factor for the consistency term calculated from pairwise
        # alignments. Valid when either of --blobalpair, --localpair,
        # --genafpair, --fastapair or --blastpair is selected. Default: 2.7
        '--weighti': ValuedParameter(Prefix='--',
                                     Name='weighti',
                                     Delimiter=' '),

        # Guide tree is built number times in the progressive stage.
        # Valid with 6mer distance. Default: 2
        '--retree': ValuedParameter(Prefix='--', Name='retree', Delimiter=' '),

        # number cycles of iterative refinement are performed. Default: 0
        '--maxiterate': ValuedParameter(Prefix='--',
                                        Name='maxiterate',
                                        Delimiter=' '),

        # Use FFT approximation in group-to-group alignment. Default: on
        '--fft': FlagParameter(Prefix='--', Name='fft'),

        # Do not use FFT approximation in group-to-group alignment.
        # Default: off
        '--nofft': FlagParameter(Prefix='--', Name='nofft'),

        # Alignment score is not checked in the iterative refinement stage.
        # Default: off (score is checked)
        '--noscore': FlagParameter(Prefix='--', Name='noscore'),

        # Use the Myers-Miller (1988) algorithm.
        # Default: automatically turned on
        # when the alignment length exceeds 10,000 (aa/nt).
        '--memsave': FlagParameter(Prefix='--', Name='memsave'),

        # Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the
        # 6mer distance. Recommended for a large number (> ~10,000) of sequences are
        # input. Default: off
        '--parttree': FlagParameter(Prefix='--', Name='parttree'),

        # The PartTree algorithm is used with distances based on DP. Slightly more
        # accurate and slower than --parttree. Recommended for a large number
        # (> ~10,000) of sequences are input. Default: off
        '--dpparttree': FlagParameter(Prefix='--', Name='dpparttree'),

        # The PartTree algorithm is used with distances based on FASTA. Slightly
        # more accurate and slower than --parttree. Recommended for a large number
        # (> ~10,000) of sequences are input. FASTA is required. Default: off
        '--fastaparttree': FlagParameter(Prefix='--', Name='fastaparttree'),

        # The number of partitions in the PartTree algorithm. Default: 50
        '--partsize': ValuedParameter(Prefix='--',
                                      Name='partsize',
                                      Delimiter=' '),

        # Do not make alignment larger than number sequences. Valid only with the
        # --*parttree options. Default: the number of input sequences
        '--groupsize': ValuedParameter(Prefix='--', Name='groupsize', Delimiter=' '),

        # Parameter

        # Gap opening penalty at group-to-group alignment. Default: 1.53
        '--op': ValuedParameter(Prefix='--', Name='op', Delimiter=' '),

        # Offset value, which works like gap extension penalty, for group-to-group
        # alignment. Deafult: 0.123
        '--ep': ValuedParameter(Prefix='--', Name='ep', Delimiter=' '),

        # Gap opening penalty at local pairwise alignment. Valid when the
        # --localpair or --genafpair option is selected. Default: -2.00
        '--lop': ValuedParameter(Prefix='--', Name='lop', Delimiter=' '),

        # Offset value at local pairwise alignment. Valid when the --localpair
        # or --genafpair option is selected. Default: 0.1
        '--lep': ValuedParameter(Prefix='--', Name='lep', Delimiter=' '),

        # Gap extension penalty at local pairwise alignment. Valid when the
        # --localpair or --genafpair option is selected. Default: -0.1
        '--lexp': ValuedParameter(Prefix='--', Name='lexp', Delimiter=' '),

        # Gap opening penalty to skip the alignment. Valid when the --genafpair
        # option is selected. Default: -6.00
        '--LOP': ValuedParameter(Prefix='--', Name='LOP', Delimiter=' '),

        # Gap extension penalty to skip the alignment. Valid when the --genafpair
        # option is selected. Default: 0.00
        '--LEXP': ValuedParameter(Prefix='--', Name='LEXP', Delimiter=' '),

        # BLOSUM number matrix (Henikoff and Henikoff 1992) is used. number=30, 45,
        # 62 or 80. Default: 62
        '--bl': ValuedParameter(Prefix='--', Name='bl', Delimiter=' '),

        # JTT PAM number (Jones et al. 1992) matrix is used. number>0.
        # Default: BLOSUM62
        '--jtt': ValuedParameter(Prefix='--', Name='jtt', Delimiter=' '),

        # Transmembrane PAM number (Jones et al. 1994) matrix is used. number>0.
        # Default: BLOSUM62
        '--tm': ValuedParameter(Prefix='--', Name='tm', Delimiter=' '),

        # Use a user-defined AA scoring matrix. The format of matrixfile is the same
        # to that of BLAST. Ignored when nucleotide sequences are input.
        # Default: BLOSUM62
        '--aamatrix': ValuedParameter(Prefix='--', Name='aamatrix', Delimiter=' '),

        # Incorporate the AA/nuc composition information into the scoring matrix.
        # Deafult: off
        '--fmodel': FlagParameter(Prefix='--', Name='fmodel'),

        # Output

        # Output format: clustal format. Default: off (fasta format)
        '--clustalout': FlagParameter(Prefix='--', Name='clustalout'),

        # Output order: same as input. Default: on
        '--inputorder': FlagParameter(Prefix='--', Name='inputorder'),

        # Output order: aligned. Default: off (inputorder)
        '--reorder': FlagParameter(Prefix='--', Name='reorder'),

        # Guide tree is output to the input.tree file. Default: off
        '--treeout': FlagParameter(Prefix='--', Name='treeout'),

        # Do not report progress. Default: off
        '--quiet': FlagParameter(Prefix='--', Name='quiet'),

        # Input

        # Assume the sequences are nucleotide. Deafult: auto
        '--nuc': FlagParameter(Prefix='--', Name='nuc'),

        # Assume the sequences are amino acid. Deafult: auto
        '--amino': FlagParameter(Prefix='--', Name='amino'),

        # Seed alignments given in alignment_n (fasta format) are aligned with
        # sequences in input. The alignment within every seed is preserved.
        '--seed': ValuedParameter(Prefix='--', Name='seed', Delimiter=' '),
    }

    _parameters = {}
    _parameters.update(_options)
    _command = "mafft"
    _suppress_stderr = False

    def _input_as_seqs(self, data):
        """Format a list of seq as input.

        Parameters
        ----------
        data: list of strings
            Each string is a sequence to be aligned.

        Returns
        -------
        A temp file name that contains the sequences.

        See Also
        --------
        burrito.util.CommandLineApplication
        """
        lines = []
        for i, s in enumerate(data):
            # will number the sequences 1,2,3,etc.
            lines.append(''.join(['>', str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)

    def _tree_out_filename(self):
        if self.Parameters['--treeout'].isOn():
            tree_filename = self._absolute(str(self._input_filename))+'.tree'
        else:
            raise ValueError("No tree output file specified.")
        return tree_filename

    def getHelp(self):
        """Method that points to the Mafft documentation."""

        help_str = ("See Mafft documentation at:\n"
                    "http://mafft.cbrc.jp/alignment/software/"
                    "mafft/software/manual/manual.html")
        return help_str

    def _get_result_paths(self, data):
        result = {}
        if self.Parameters['--treeout'].isOn():
            out_name = self._tree_out_filename()
            result['Tree'] = ResultPath(Path=out_name, IsWritten=True)
        return result


def align_unaligned_seqs(seqs_fp, moltype=DNA, params=None, accurate=False):
    """Aligns unaligned sequences

    Parameters
    ----------
    seqs_fp : string
        file path of the input fasta file
    moltype : {skbio.DNA, skbio.RNA, skbio.Protein}
    params : dict-like type
        It pass the additional parameter settings to the application.
        Default is None.
    accurate : boolean
        Perform accurate alignment or not. It will sacrifice performance
        if set to True. Default is False.

    Returns
    -------
    Alignment object
        The aligned sequences.

    See Also
    --------
    skbio.Alignment
    skbio.DNA
    skbio.RNA
    skbio.Protein
    """
    # Create Mafft app.
    app = Mafft(InputHandler='_input_as_path', params=params)

    # Turn on correct sequence type
    app.Parameters[MOLTYPE_MAP[moltype]].on()

    # Do not report progress
    app.Parameters['--quiet'].on()

    # More accurate alignment, sacrificing performance.
    if accurate:
        app.Parameters['--globalpair'].on()
        app.Parameters['--maxiterate'].Value = 1000

    # Get results using int_map as input to app
    res = app(seqs_fp)

    # Get alignment as dict out of results
    alignment = Alignment.read(res['StdOut'], constructor=moltype)

    # Clean up
    res.cleanUp()

    return alignment


def add_seqs_to_alignment(seqs_fp, aln_fp, moltype=DNA,
                          params=None, accurate=False):
    """Returns an Alignment object from seqs and existing Alignment.

    The "--seed" option can be used for adding unaligned sequences into
    a highly reliable alignment (seed) consisting of a small number of
    sequences.

    Parameters
    ----------
    seqs_fp : string
        file path of the unaligned sequences
    aln_fp : string
        file path of the seed alignment
    params : dict of parameters to pass in to the Mafft app controller.

    Returns
    -------
        The aligned sequences. The seq in the seed alignment will have
        "_seed_" prefixed to their seq id.
    """
    if params is None:
        params = {'--seed': aln_fp}
    else:
        params['--seed'] = aln_fp

    return align_unaligned_seqs(seqs_fp, moltype, params, accurate)


def align_two_alignments(aln1_fp, aln2_fp, moltype, params=None):
    """Returns an Alignment object from two existing Alignments.

    Parameters
    ----------
    aln1_fp : string
        file path of 1st alignment
    aln2_fp : string
        file path of 2nd alignment
    params : dict of parameters to pass in to the Mafft app controller.

    Returns
    -------
        The aligned sequences.
    """

    # Create Mafft app.
    app = Mafft(InputHandler='_input_as_paths',
                params=params,
                SuppressStderr=False)
    app._command = 'mafft-profile'

    # Get results using int_map as input to app
    res = app([aln1_fp, aln2_fp])

    return Alignment.read(res['StdOut'], constructor=moltype)
