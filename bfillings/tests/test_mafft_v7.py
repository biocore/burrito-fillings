#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from os import getcwd
from os.path import join
import tempfile
import shutil
from unittest import TestCase, main
from skbio import RNA
from bfillings.mafft_v7 import (Mafft, align_unaligned_seqs,
                                add_seqs_to_alignment, align_two_alignments)


class MafftTests(TestCase):
    """Tests for the Mafft application controller"""
    def setUp(self):
        """Mafft general setUp method for all tests"""
        self.temp_dir = tempfile.mkdtemp()

        self.seqs1 = ('>1\n'
                      'ACUGCUAGCUAGUAGCGUACGUA\n'
                      '>2\n'
                      'GCUACGUAGCUAC\n'
                      '>3\n'
                      'GCGGCUAUUAGAUCGUA\n')
        self.seqs1_fp = join(self.temp_dir, 'seq1.fa')
        with open(self.seqs1_fp, 'w') as f:
            f.write(self.seqs1)
        self.seqs1_aln = ('>1\n---acugcuagcuaguagcguacgua\n'
                          '>2\n------gcuacguagcuac-------\n'
                          '>3\ngcggcuauuagaucgua---------\n')
        self.seqs1_aln_fp = join(self.temp_dir, 'seq1_aln.fa')
        with open(self.seqs1_aln_fp, 'w') as f:
            f.write(self.seqs1_aln)

        self.seqs2 = ('>a\nUAGGCUCUGAUAUAAUAGCUCUC\n'
                      '>b\nUAUCGCUUCGACGAUUCUCUGAUAGAGA\n'
                      '>c\nUGACUACGCAU\n')
        self.seqs2_fp = join(self.temp_dir, 'seq2.fa')
        with open(self.seqs2_fp, 'w') as f:
            f.write(self.seqs2)

        self.add_seqs_aligned = (">_seed_1\n"
                                 "----------acugcuagcuaguagcguacgua\n"
                                 ">_seed_2\n"
                                 "-------------gcuacguagcuac-------\n"
                                 ">_seed_3\n"
                                 "-------gcggcuauuagaucgua---------\n"
                                 ">a\n"
                                 "-------uaggcucugauauaauagcucuc---\n"
                                 ">b\n"
                                 "uaucgcuucgacgauucucugauagaga-----\n"
                                 ">c\n"
                                 "-------------------ugacuacgcau---\n")

        self.align1 = (">seq_0\nACUGCUAGCUAGUAGCGUACGUA\n"
                       ">seq_1\nGCUACGUAGCUAC----------\n"
                       ">seq_2\nGCGGCUAUUAGAU------CGUA\n")
        self.align1_fp = join(self.temp_dir, 'align1.fa')
        with open(self.align1_fp, 'w') as f:
            f.write(self.align1)
        self.align2 = (">a\nUAGGCUCUGAUAUAAUAGCUCUC---------\n"
                       ">b\nUA----UCGCUUCGACGAUUCUCUGAUAGAGA\n"
                       ">c\nUG------------ACUACGCAU---------\n")
        self.align2_fp = join(self.temp_dir, 'align2.fa')
        with open(self.align2_fp, 'w') as f:
            f.write(self.align2)
        self.align_two_align = (">seq_0\n"
                                "--------------acugcuagcuaguagcguacgua\n"
                                ">seq_1\n"
                                "--------------gcuacguagcuac----------\n"
                                ">seq_2\n"
                                "--------------gcggcuauuagau------cgua\n"
                                ">a\n"
                                "uaggcucugauauaauagcucuc--------------\n"
                                ">b\n"
                                "ua----ucgcuucgacgauucucugauagaga-----\n"
                                ">c\n"
                                "ug------------acuacgcau--------------\n")

    def test_base_command(self):
        """Mafft BaseCommand should return the correct BaseCommand"""
        c = Mafft()
        self.assertEqual(c.BaseCommand,
                         ''.join(['cd "', getcwd(), '/"; ', 'mafft']))
        c.Parameters['--quiet'].on()
        self.assertEqual(c.BaseCommand,
                         ''.join(['cd "', getcwd(), '/"; ', 'mafft --quiet']))
        c.Parameters['--globalpair'].on()
        self.assertEqual(c.BaseCommand,
                         'cd "%s/"; mafft --globalpair --quiet' % getcwd())
        c.Parameters['--maxiterate'].on(1000)
        self.assertEqual(
            c.BaseCommand,
            ('cd "%s/"; '
             'mafft --maxiterate 1000 --globalpair --quiet' % getcwd()))

    def test_align_unaligned_seqs(self):
        """align_unaligned_seqs should work as expected"""
        res = align_unaligned_seqs(self.seqs1_fp, RNA)
        self.assertEqual(res.toFasta(), self.seqs1_aln)

    def test_add_seqs_to_alignment(self):
        """add_seqs_to_alignment should work as expected."""
        res = add_seqs_to_alignment(self.seqs2_fp, self.seqs1_aln_fp, RNA)
        self.assertEqual(res.toFasta(), self.add_seqs_aligned)

    def test_align_two_alignments(self):
        """align_two_alignments should work as expected."""
        res = align_two_alignments(self.align1_fp, self.align2_fp, RNA)
        self.assertEqual(res.toFasta(), self.align_two_align)

    def tearDown(self):
        """Last test executed: cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)


if __name__ == '__main__':
    main()
