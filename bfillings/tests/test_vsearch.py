#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2015--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
"""
Unit tests for the VSEARCH version 1.1.1 Application controller
===============================================================
"""


from unittest import TestCase, main
from os import close
from os.path import exists, join, dirname
from tempfile import mkstemp, mkdtemp
from shutil import rmtree

from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta

from bfillings.vsearch import (vsearch_dereplicate_exact_seqs,
                               vsearch_sort_by_abundance,
                               vsearch_chimera_filter_de_novo,
                               vsearch_chimera_filter_ref)


# Test class and cases
class VsearchTests(TestCase):
    """ Tests for VSEARCH version 1.1.1 functionality """

    def setUp(self):
        self.output_dir = mkdtemp()
        self.seqs_to_derep = seqs_to_derep
        self.seqs_to_derep_max_min_abundance =\
            seqs_to_derep_max_min_abundance
        self.seqs_to_derep_merged_derep_files =\
            seqs_to_derep_merged_derep_files
        self.seqs_to_sort = seqs_to_sort
        self.amplicon_reads = amplicon_reads
        self.single_chimera = single_chimera
        self.single_chimera_ref = single_chimera_ref
        self.uchime_ref_db = uchime_ref_db
        self.uchime_single_ref_db = uchime_single_ref_db
        self.uc_output = uc_output
        self.uc_output_incorrect_1 = uc_output_incorrect_1

        # temporary file for seqs_to_derep
        f, self.seqs_to_derep_fp = mkstemp(prefix='tmp_seqs_to_derep_',
                                           suffix='.fasta')
        close(f)
        # write seqs_to_derep to file
        with open(self.seqs_to_derep_fp, 'w') as tmp:
            tmp.write(self.seqs_to_derep)

        # temporary file for seqs_to_derep_max_min_abundance
        f, self.seqs_to_derep_max_min_abundance_fp =\
            mkstemp(prefix='tmp_seqs_to_derep_abun_',
                    suffix='.fasta')
        close(f)
        # write seqs_to_derep_max_min_abundance to file
        with open(self.seqs_to_derep_max_min_abundance_fp, 'w') as tmp:
            tmp.write(self.seqs_to_derep_max_min_abundance)

        # temporary file for seqs_to_derep_merged_derep_files
        f, self.seqs_to_derep_merged_derep_files_fp =\
            mkstemp(prefix='tmp_seqs_to_derep_concat_',
                    suffix='.fasta')
        close(f)
        # write seqs_to_derep_merged_derep_files to file
        with open(self.seqs_to_derep_merged_derep_files_fp, 'w') as tmp:
            tmp.write(self.seqs_to_derep_merged_derep_files)

        # temporary file for seqs_to_sort
        f, self.seqs_to_sort_fp = mkstemp(prefix='tmp_seqs_to_sort_',
                                          suffix='.fasta')
        close(f)
        # write seqs_to_sort to file
        with open(self.seqs_to_sort_fp, 'w') as tmp:
            tmp.write(self.seqs_to_sort)

        # temporary file for amplicon_reads
        f, self.amplicon_reads_fp = mkstemp(prefix='tmp_amplicon_reads_',
                                            suffix='.fasta')
        close(f)
        # write amplicon_reads to file
        with open(self.amplicon_reads_fp, 'w') as tmp:
            tmp.write(self.amplicon_reads)

        # temporary file for single_chimera
        f, self.single_chimera_fp = mkstemp(prefix='tmp_single_chimera_',
                                            suffix='.fasta')
        close(f)
        # write single_chimera to file
        # (de novo chimera checking)
        with open(self.single_chimera_fp, 'w') as tmp:
            tmp.write(self.single_chimera)

        # temporary file for single_chimera_ref
        f, self.single_chimera_ref_fp = mkstemp(prefix='tmp_single_chimera_',
                                                suffix='.fasta')
        close(f)
        # write single_chimera_ref to file
        # (reference chimera checking)
        with open(self.single_chimera_ref_fp, 'w') as tmp:
            tmp.write(self.single_chimera_ref)

        # temporary file for uchime_ref_db
        f, self.uchime_ref_db_fp = mkstemp(prefix='tmp_uchime_ref_db_',
                                           suffix='.fasta')
        close(f)
        # write uchime_ref_db to file
        with open(self.uchime_ref_db_fp, 'w') as tmp:
            tmp.write(self.uchime_ref_db)

        # temporary file for uchime_single_ref_db
        f, self.uchime_single_ref_db_fp =\
            mkstemp(prefix='tmp_uchime_single_ref_db_',
                    suffix='.fasta')
        close(f)
        # write uchime_single_ref_db to file
        with open(self.uchime_single_ref_db_fp, 'w') as tmp:
            tmp.write(self.uchime_single_ref_db)

        # temporary file for .uc output
        f, self.uc_fp =\
            mkstemp(prefix='tmp_uc_',
                    suffix='.uc')
        close(f)
        # write uc_output to file
        with open(self.uc_fp, 'w') as tmp:
            tmp.write(self.uc_output)

        # temporary file for .uc incorrect output 1
        f, self.uc_incorrect_1_fp =\
            mkstemp(prefix='tmp_uc_incorrect_1',
                    suffix='.uc')
        close(f)
        # write uc_output_incorrect_1 to file
        with open(self.uc_incorrect_1_fp, 'w') as tmp:
            tmp.write(self.uc_output_incorrect_1)

        # list of files to remove
        self.files_to_remove = [self.seqs_to_derep_fp,
                                self.seqs_to_derep_max_min_abundance_fp,
                                self.seqs_to_derep_merged_derep_files_fp,
                                self.seqs_to_sort_fp,
                                self.amplicon_reads_fp,
                                self.single_chimera_fp,
                                self.single_chimera_ref_fp,
                                self.uchime_ref_db_fp,
                                self.uchime_single_ref_db_fp,
                                self.uc_fp,
                                self.uc_incorrect_1_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.output_dir)

    def test_vsearch_chimera_filter_ref(self):
        """ Test reference chimera filter, output only
            chimeric sequences and log
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_ref(
                self.amplicon_reads_fp,
                self.output_dir,
                self.uchime_ref_db_fp,
                output_chimeras=True,
                output_nonchimeras=False,
                output_alns=False,
                output_tabular=False,
                log_name="vsearch_uchime_ref_chimera_filtering.log",
                HALT_EXEC=False)

        self.assertTrue(nonchimeras_fp is None)
        self.assertTrue(alns_fp is None)
        self.assertTrue(tabular_fp is None)
        self.assertTrue(exists(log_fp))

        expected_chimeras = ['251;size=2;', '320;size=2;', '36;size=2;',
                             '672;size=2;', '142;size=1;', '201;size=1;',
                             '241;size=1;', '279;size=1;', '299;size=1;',
                             '359;size=1;', '375;size=1;', '407;size=1;',
                             '423;size=1;', '516;size=1;', '618;size=1;',
                             '717;size=1;', '902;size=1;', '918;size=1;',
                             '941;size=1;']

        num_seqs = 0

        with open(chimeras_fp, "U") as chimeras_f:
            for label, seq in parse_fasta(chimeras_f):
                # check label represents chimeric sequence
                self.assertTrue(label in expected_chimeras)
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 19)

    def test_vsearch_chimera_filter_ref_output(self):
        """ Raise error when no output is selected for
            reference chimera filtering
        """

        self.assertRaises(ValueError,
                          vsearch_chimera_filter_ref,
                          fasta_filepath=self.amplicon_reads_fp,
                          working_dir=self.output_dir,
                          db_filepath=self.uchime_ref_db_fp,
                          output_chimeras=False,
                          output_nonchimeras=False,
                          output_alns=False,
                          output_tabular=False,
                          log_name="vsearch_uchime_ref_chimera_filtering.log",
                          HALT_EXEC=False)

    def test_vsearch_chimera_filter_ref_output_nonchimeras(self):
        """ Test ref chimera filter, output nonchimeric sequences
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_ref(
                self.amplicon_reads_fp,
                self.output_dir,
                self.uchime_ref_db_fp,
                output_chimeras=False,
                output_nonchimeras=True,
                output_alns=False,
                output_tabular=False,
                log_name="vsearch_uchime_ref_chimera_filtering.log",
                HALT_EXEC=False)

        self.assertTrue(chimeras_fp is None)
        self.assertTrue(alns_fp is None)
        self.assertTrue(tabular_fp is None)
        self.assertTrue(exists(log_fp))

        expected_nonchimeras =\
            ['3;size=102;', '16;size=95;', '22;size=93;', '2;size=87;',
             '39;size=84;', '4;size=79;', '6;size=72;', '11;size=70;',
             '45;size=67;', '1;size=65;', '425;size=2;', '100;size=1;',
             '102;size=1;', '10;size=1;', '115;size=1;', '123;size=1;',
             '132;size=1;', '134;size=1;', '140;size=1;', '144;size=1;',
             '148;size=1;', '14;size=1;', '156;size=1;', '15;size=1;',
             '161;size=1;', '162;size=1;', '186;size=1;', '203;size=1;',
             '217;size=1;', '218;size=1;', '21;size=1;', '221;size=1;',
             '222;size=1;', '225;size=1;', '233;size=1;', '234;size=1;',
             '235;size=1;', '249;size=1;', '24;size=1;', '259;size=1;',
             '266;size=1;', '26;size=1;', '27;size=1;', '296;size=1;',
             '303;size=1;', '306;size=1;', '307;size=1;', '322;size=1;',
             '326;size=1;', '32;size=1;', '332;size=1;', '333;size=1;',
             '338;size=1;', '360;size=1;', '362;size=1;', '364;size=1;',
             '366;size=1;', '369;size=1;', '371;size=1;', '373;size=1;',
             '374;size=1;', '37;size=1;', '386;size=1;', '387;size=1;',
             '392;size=1;', '393;size=1;', '397;size=1;', '405;size=1;',
             '414;size=1;', '418;size=1;', '431;size=1;', '436;size=1;',
             '444;size=1;', '445;size=1;', '456;size=1;', '460;size=1;',
             '469;size=1;', '470;size=1;', '477;size=1;', '479;size=1;',
             '486;size=1;', '500;size=1;', '515;size=1;', '528;size=1;',
             '530;size=1;', '531;size=1;', '549;size=1;', '551;size=1;',
             '557;size=1;', '559;size=1;', '561;size=1;', '562;size=1;',
             '564;size=1;', '566;size=1;', '568;size=1;', '570;size=1;',
             '578;size=1;', '57;size=1;', '586;size=1;', '596;size=1;',
             '600;size=1;', '612;size=1;', '625;size=1;', '632;size=1;',
             '649;size=1;', '650;size=1;', '651;size=1;', '664;size=1;',
             '66;size=1;', '673;size=1;', '675;size=1;', '682;size=1;',
             '690;size=1;', '699;size=1;', '709;size=1;', '73;size=1;',
             '740;size=1;', '745;size=1;', '746;size=1;', '748;size=1;',
             '760;size=1;', '766;size=1;', '778;size=1;', '77;size=1;',
             '791;size=1;', '797;size=1;', '7;size=1;', '809;size=1;',
             '813;size=1;', '814;size=1;', '816;size=1;', '817;size=1;',
             '821;size=1;', '824;size=1;', '827;size=1;', '82;size=1;',
             '83;size=1;', '842;size=1;', '851;size=1;', '853;size=1;',
             '862;size=1;', '863;size=1;', '866;size=1;', '871;size=1;',
             '879;size=1;', '886;size=1;', '892;size=1;', '895;size=1;',
             '897;size=1;', '904;size=1;', '912;size=1;', '916;size=1;',
             '91;size=1;', '920;size=1;', '921;size=1;', '925;size=1;',
             '930;size=1;', '942;size=1;', '945;size=1;', '947;size=1;',
             '948;size=1;', '952;size=1;', '956;size=1;', '958;size=1;',
             '964;size=1;', '967;size=1;', '984;size=1;', '992;size=1;',
             '993;size=1;']

        num_seqs = 0

        # check nonchimeras fasta file
        with open(nonchimeras_fp, "U") as nonchimeras_f:
            for label, seq in parse_fasta(nonchimeras_f):
                # check label represents chimeric sequence
                self.assertTrue(label in expected_nonchimeras)
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 169)

    def test_vsearch_chimera_filter_ref_output_alns_tab(self):
        """ Test ref chimera filter, output only
            chimeric alignments and tabular format
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_ref(
                self.single_chimera_ref_fp,
                self.output_dir,
                self.uchime_single_ref_db_fp,
                output_chimeras=False,
                output_nonchimeras=False,
                output_alns=True,
                output_tabular=True,
                log_name="vsearch_uchime_ref_chimera_filtering.log",
                HALT_EXEC=False)

        self.assertTrue(chimeras_fp is None)
        self.assertTrue(nonchimeras_fp is None)
        self.assertTrue(exists(log_fp))

        # check alignment is correct
        with open(alns_fp, 'U') as alns_f:
            actual_alns = alns_f.read()
        self.assertEquals(single_chimera_ref_aln, actual_alns)

        # check tabular output is correct
        with open(tabular_fp, 'U') as tabular_f:
            actual_tab = tabular_f.read()

        self.assertEquals(single_chimera_ref_tab, actual_tab)

    def test_vsearch_chimera_filter_de_novo(self):
        """ Test de novo chimera filter, output only
            chimeric sequences
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_de_novo(
                self.amplicon_reads_fp,
                self.output_dir,
                output_chimeras=True,
                output_nonchimeras=False,
                output_alns=False,
                output_tabular=False,
                log_name="vsearch_uchime_de_novo_chimera_filtering.log",
                HALT_EXEC=False)

        self.assertTrue(nonchimeras_fp is None)
        self.assertTrue(alns_fp is None)
        self.assertTrue(tabular_fp is None)
        self.assertTrue(exists(log_fp))

        expected_chimeras = ['251;size=2;', '320;size=2;', '36;size=2;',
                             '672;size=2;', '142;size=1;', '201;size=1;',
                             '241;size=1;', '279;size=1;', '299;size=1;',
                             '359;size=1;', '375;size=1;', '407;size=1;',
                             '423;size=1;', '516;size=1;', '618;size=1;',
                             '717;size=1;', '902;size=1;', '918;size=1;',
                             '941;size=1;']

        num_seqs = 0

        with open(chimeras_fp, "U") as chimeras_f:
            for label, seq in parse_fasta(chimeras_f):
                # check label represents chimeric sequence
                self.assertTrue(label in expected_chimeras)
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 19)

    def test_vsearch_chimera_filter_de_novo_output(self):
        """ Raise error when no output is selected for
            de novo chimera filtering
        """

        self.assertRaises(ValueError,
                          vsearch_chimera_filter_de_novo,
                          fasta_filepath=self.amplicon_reads_fp,
                          working_dir=self.output_dir,
                          output_chimeras=False,
                          output_nonchimeras=False,
                          output_alns=False,
                          output_tabular=False,
                          log_name="vsearch_uchime_de_novo_chimera_filter.log",
                          HALT_EXEC=False)

    def test_vsearch_chimera_filter_de_novo_output_nonchimeras(self):
        """ Test de novo chimera filter, output nonchimeric sequences
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_de_novo(
                self.amplicon_reads_fp,
                self.output_dir,
                output_chimeras=False,
                output_nonchimeras=True,
                output_alns=False,
                output_tabular=False,
                log_name="vsearch_uchime_de_novo_chimera_filter.log",
                HALT_EXEC=False)

        self.assertTrue(chimeras_fp is None)
        self.assertTrue(alns_fp is None)
        self.assertTrue(tabular_fp is None)
        self.assertTrue(exists(log_fp))

        expected_nonchimeras =\
            ['3;size=102;', '16;size=95;', '22;size=93;', '2;size=87;',
             '39;size=84;', '4;size=79;', '6;size=72;', '11;size=70;',
             '45;size=67;', '1;size=65;', '425;size=2;', '100;size=1;',
             '102;size=1;', '10;size=1;', '115;size=1;', '123;size=1;',
             '132;size=1;', '134;size=1;', '140;size=1;', '144;size=1;',
             '148;size=1;', '14;size=1;', '156;size=1;', '15;size=1;',
             '161;size=1;', '162;size=1;', '186;size=1;', '203;size=1;',
             '217;size=1;', '218;size=1;', '21;size=1;', '221;size=1;',
             '222;size=1;', '225;size=1;', '233;size=1;', '234;size=1;',
             '235;size=1;', '249;size=1;', '24;size=1;', '259;size=1;',
             '266;size=1;', '26;size=1;', '27;size=1;', '296;size=1;',
             '303;size=1;', '306;size=1;', '307;size=1;', '322;size=1;',
             '326;size=1;', '32;size=1;', '332;size=1;', '333;size=1;',
             '338;size=1;', '360;size=1;', '362;size=1;', '364;size=1;',
             '366;size=1;', '369;size=1;', '371;size=1;', '373;size=1;',
             '374;size=1;', '37;size=1;', '386;size=1;', '387;size=1;',
             '392;size=1;', '393;size=1;', '397;size=1;', '405;size=1;',
             '414;size=1;', '418;size=1;', '431;size=1;', '436;size=1;',
             '444;size=1;', '445;size=1;', '456;size=1;', '460;size=1;',
             '469;size=1;', '470;size=1;', '477;size=1;', '479;size=1;',
             '486;size=1;', '500;size=1;', '515;size=1;', '528;size=1;',
             '530;size=1;', '531;size=1;', '549;size=1;', '551;size=1;',
             '557;size=1;', '559;size=1;', '561;size=1;', '562;size=1;',
             '564;size=1;', '566;size=1;', '568;size=1;', '570;size=1;',
             '578;size=1;', '57;size=1;', '586;size=1;', '596;size=1;',
             '600;size=1;', '612;size=1;', '625;size=1;', '632;size=1;',
             '649;size=1;', '650;size=1;', '651;size=1;', '664;size=1;',
             '66;size=1;', '673;size=1;', '675;size=1;', '682;size=1;',
             '690;size=1;', '699;size=1;', '709;size=1;', '73;size=1;',
             '740;size=1;', '745;size=1;', '746;size=1;', '748;size=1;',
             '760;size=1;', '766;size=1;', '778;size=1;', '77;size=1;',
             '791;size=1;', '797;size=1;', '7;size=1;', '809;size=1;',
             '813;size=1;', '814;size=1;', '816;size=1;', '817;size=1;',
             '821;size=1;', '824;size=1;', '827;size=1;', '82;size=1;',
             '83;size=1;', '842;size=1;', '851;size=1;', '853;size=1;',
             '862;size=1;', '863;size=1;', '866;size=1;', '871;size=1;',
             '879;size=1;', '886;size=1;', '892;size=1;', '895;size=1;',
             '897;size=1;', '904;size=1;', '912;size=1;', '916;size=1;',
             '91;size=1;', '920;size=1;', '921;size=1;', '925;size=1;',
             '930;size=1;', '942;size=1;', '945;size=1;', '947;size=1;',
             '948;size=1;', '952;size=1;', '956;size=1;', '958;size=1;',
             '964;size=1;', '967;size=1;', '984;size=1;', '992;size=1;',
             '993;size=1;']

        num_seqs = 0

        # check nonchimeras fasta file
        with open(nonchimeras_fp, "U") as nonchimeras_f:
            for label, seq in parse_fasta(nonchimeras_f):
                # check label represents chimeric sequence
                self.assertTrue(label in expected_nonchimeras)
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 169)

    def test_vsearch_chimera_filter_de_novo_output_alns_tab(self):
        """ Test de novo chimera filter, output only
            chimeric alignments and tabular format
        """
        chimeras_fp, nonchimeras_fp, alns_fp, tabular_fp, log_fp =\
            vsearch_chimera_filter_de_novo(
                self.single_chimera_fp,
                self.output_dir,
                output_chimeras=False,
                output_nonchimeras=False,
                output_alns=True,
                output_tabular=True,
                log_name="vsearch_uchime_de_novo_chimera_filter.log",
                HALT_EXEC=False)

        self.assertTrue(chimeras_fp is None)
        self.assertTrue(nonchimeras_fp is None)
        self.assertTrue(exists(log_fp))

        # check alignment is correct
        with open(alns_fp, 'U') as alns_f:
            actual_alns = alns_f.read()
        self.assertEquals(single_chimera_aln, actual_alns)

        # check tabular output is correct
        with open(tabular_fp, 'U') as tabular_f:
            actual_tab = tabular_f.read()
        self.assertEquals(single_chimera_tab, actual_tab)

    def test_vsearch_sort_by_abundance(self):
        """ Test sorting sequences by abundance
        """
        tmp_fp = join(self.output_dir, "tmp_sorted_reads.fasta")

        output_sorted, log_fp = vsearch_sort_by_abundance(
            self.seqs_to_sort_fp,
            tmp_fp,
            working_dir=None,
            minsize=None,
            maxsize=None,
            log_name="abundance_sort.log",
            HALT_EXEC=False)

        self.assertTrue(exists(log_fp))

        expected_order = ['HWI-ST157_0368:1:2107:19923:3944#0/1;size=100;',
                          'HWI-ST157_0368:1:1201:8401:113582#0/1;size=10;',
                          'HWI-ST157_0368:1:2204:20491:181552#0/1;size=10;',
                          'HWI-ST157_0368:1:2105:3428:36721#0/1;size=5;',
                          'HWI-ST157_0368:1:2105:6731:137157#0/1;size=4;',
                          'HWI-ST157_0368:1:2106:18272:88408#0/1;size=2;',
                          'HWI-ST157_0368:1:1106:12200:200135#0/1;size=1;',
                          'HWI-ST157_0368:1:2208:9135:145970#0/1;size=1;']

        num_seqs = 0

        with open(output_sorted, "U") as tmp_f:
            for label, seq in parse_fasta(tmp_f):
                # check label is in correct order
                self.assertEquals(label, expected_order[num_seqs])
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 8)

    def test_vsearch_sort_by_abundance_minsize_1_maxsize_10(self):
        """ Test sorting sequences by abundance,
            discard sequences with an abundance value smaller
            than 1 and greater than 10
        """
        tmp_fp = join(self.output_dir, "tmp_sorted_reads.fasta")

        output_sorted, log_fp = vsearch_sort_by_abundance(
            self.seqs_to_sort_fp,
            tmp_fp,
            working_dir=None,
            minsize=2,
            maxsize=10,
            log_name="abundance_sort.log",
            HALT_EXEC=False)

        self.assertTrue(exists(log_fp))

        expected_order = ['HWI-ST157_0368:1:1201:8401:113582#0/1;size=10;',
                          'HWI-ST157_0368:1:2204:20491:181552#0/1;size=10;',
                          'HWI-ST157_0368:1:2105:3428:36721#0/1;size=5;',
                          'HWI-ST157_0368:1:2105:6731:137157#0/1;size=4;',
                          'HWI-ST157_0368:1:2106:18272:88408#0/1;size=2;']

        num_seqs = 0

        with open(output_sorted, "U") as tmp_f:
            for label, seq in parse_fasta(tmp_f):
                # check label is in correct order
                self.assertEquals(label, expected_order[num_seqs])
                # check sequence exists
                self.assertTrue(len(seq) > 0)
                num_seqs += 1

        self.assertTrue(num_seqs, 5)

    def test_vsearch_dereplicate_exact_seqs(self):
    	""" Test dereplicating sequences
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        dereplicated_seqs_fp, uc_fp, log_fp = vsearch_dereplicate_exact_seqs(
            self.seqs_to_derep_fp,
            tmp_fp,
            output_uc=False,
            working_dir=self.output_dir,
            strand="both",
            maxuniquesize=None,
            minuniquesize=None,
            sizein=False,
            sizeout=True)

        # no output for .uc
        self.assertTrue(uc_fp is None)
        self.assertTrue(exists(log_fp))

        num_seqs = 0
        expected_derep = ['HWI-ST157_0368:1:1207:16180:126921#0/1;size=3;',
                          'HWI-ST157_0368:1:2103:7895:197066#0/1;size=3;',
                          'HWI-ST157_0368:1:1106:11378:83198#0/1;size=1;',
                          'HWI-ST157_0368:1:2102:15078:69955#0/1;size=1;']

        with open(tmp_fp, "U") as tmp_f:
            for label, seq in parse_fasta(tmp_f):
                num_seqs += 1
                # check output labels are correct
                self.assertTrue(label in expected_derep)
                # check sequence exists
                self.assertTrue(len(seq) > 0)

        # check there are 4 sequences after dereplication
        self.assertEquals(num_seqs, 4)

    def test_vsearch_dereplicate_exact_seqs_uc(self):
        """ Test dereplicating sequences with .uc output
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        dereplicated_seqs_fp, uc_fp, log_fp = vsearch_dereplicate_exact_seqs(
            self.seqs_to_derep_fp,
            tmp_fp,
            output_uc=True,
            working_dir=self.output_dir,
            strand="both",
            maxuniquesize=None,
            minuniquesize=None,
            sizein=False,
            sizeout=True)

        # .uc exists
        self.assertTrue(exists(uc_fp))
        self.assertTrue(exists(log_fp))

        id_to_count = {}

        num_seqs = 0
        expected_derep = {'HWI-ST157_0368:1:1207:16180:126921#0/1': 3,
                          'HWI-ST157_0368:1:2103:7895:197066#0/1': 3,
                          'HWI-ST157_0368:1:1106:11378:83198#0/1': 1,
                          'HWI-ST157_0368:1:2102:15078:69955#0/1': 1}

        with open(uc_fp, 'U') as uc_f:
            for line in uc_f:
                if line.startswith('S'):
                    num_seqs += 1
                    label = line.strip().split('\t')[8]
                    # check output labels are correct
                    self.assertTrue(label in expected_derep)
                    id_to_count[label] = 1
                elif line.startswith('H'):
                    seed = line.strip().split('\t')[9]
                    id_to_count[seed] += 1

        # check there are 4 sequences after dereplication
        self.assertEquals(num_seqs, 4)

        for label in id_to_count:
            self.assertEquals(expected_derep[label], id_to_count[label])

    def test_vsearch_dereplicate_exact_seqs_empty_working_dir(self):
        """ Test dereplicating sequences without passing
            a working directory
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        dereplicated_seqs_fp, uc_fp, log_fp = vsearch_dereplicate_exact_seqs(
            self.seqs_to_derep_fp,
            tmp_fp,
            output_uc=True,
            working_dir=None,
            strand="both",
            maxuniquesize=None,
            minuniquesize=None,
            sizein=False,
            sizeout=True)

        self.assertTrue(exists(log_fp))

        # check dereplicated seqs and uc file in the same
        # directory (same path as tmp_fp)
        self.assertEquals(dirname(tmp_fp), dirname(dereplicated_seqs_fp))
        self.assertEquals(dirname(tmp_fp), dirname(uc_fp))

    def test_vsearch_dereplicate_exact_seqs_abundance(self):
        """ Test dereplicating sequences and discarding those with
            abundance < 2 and abundance > 6
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        dereplicated_seqs_fp, uc_fp, log_fp = vsearch_dereplicate_exact_seqs(
            self.seqs_to_derep_max_min_abundance_fp,
            tmp_fp,
            output_uc=False,
            working_dir=self.output_dir,
            strand="both",
            maxuniquesize=6,
            minuniquesize=2,
            sizein=False,
            sizeout=True)

        # no output for .uc
        self.assertTrue(uc_fp is None)
        self.assertTrue(exists(log_fp))

        num_seqs = 0
        expected_derep = ['HWI-ST157_0368:1:1106:10560:153880#0/1;size=6;',
                          'HWI-ST157_0368:1:2103:12440:90119#0/1;size=2;',
                          'HWI-ST157_0368:1:1106:15269:103850#0/1;size=3;',
                          'HWI-ST157_0368:1:1205:9745:86166#0/1;size=5;']

        with open(tmp_fp, "U") as tmp_f:
            for label, seq in parse_fasta(tmp_f):
                num_seqs += 1
                # check output labels are correct
                self.assertTrue(label in expected_derep)
                # check sequence exists
                self.assertTrue(len(seq) > 0)

        # check there are 4 sequences after dereplication
        self.assertEquals(num_seqs, 4)

    def test_vsearch_dereplicate_exact_seqs_merged(self):
        """ Test dereplicating sequences which already contain
            abundance information in the label from previous
            dereplication (ex. two dereplicated files have been
            merged into a new file for dereplication)
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        dereplicated_seqs_fp, uc_fp, log_fp = vsearch_dereplicate_exact_seqs(
            self.seqs_to_derep_merged_derep_files_fp,
            tmp_fp,
            output_uc=False,
            working_dir=self.output_dir,
            strand="both",
            maxuniquesize=None,
            minuniquesize=None,
            sizein=True,
            sizeout=True)

        # no output for .uc
        self.assertTrue(uc_fp is None)
        self.assertTrue(exists(log_fp))

        num_seqs = 0
        expected_derep = ['HWI-ST157_0368:1:1207:16180:126921#0/1;size=6;',
                          'HWI-ST157_0368:1:2103:7895:197066#0/1;size=6;',
                          'HWI-ST157_0368:1:1106:11378:83198#0/1;size=2;',
                          'HWI-ST157_0368:1:2102:15078:69955#0/1;size=2;']

        with open(tmp_fp, "U") as tmp_f:
            for label, seq in parse_fasta(tmp_f):
                num_seqs += 1
                # check output labels are correct
                self.assertTrue(label in expected_derep)
                # check sequence exists
                self.assertTrue(len(seq) > 0)

        # check there are 4 sequences after dereplication
        self.assertEquals(num_seqs, 4)

    def test_vsearch_dereplicate_exact_seqs_strand(self):
        """ Raise error when strand parameter is something
            other than 'plus' or 'both'
        """
        tmp_fp = join(self.output_dir, "tmp_derep_reads.fasta")

        self.assertRaises(ValueError,
                          vsearch_dereplicate_exact_seqs,
                          fasta_filepath=self.seqs_to_derep_fp,
                          output_filepath=tmp_fp,
                          output_uc=False,
                          working_dir=None,
                          strand="minus",
                          maxuniquesize=None,
                          minuniquesize=None,
                          sizein=False,
                          sizeout=True,
                          log_name="derep.log",
                          HALT_EXEC=False)

    def test_clusters_from_uc_file_vsearch(self):
        """ Test clusters_from_uc_file() with VSEARCH output
        """
        clusters, failures, seeds = clusters_from_uc_file(self.uc_fp)
        expected_clusters = {'s1_80': ['s1_80', 's1_81', 's1_82'],
                             's1_0': ['s1_0', 's1_1'],
                             's1_10': ['s1_10', 's1_12', 's1_13'],
                             's1_118': ['s1_118', 's1_119'],
                             's1_11': ['s1_11'],
                             's1_25': ['s1_25']}
        self.assertDictEqual(clusters, expected_clusters)

    def test_clusters_from_uc_file_vsearch_duplicate_seed_ids(self):
        """ Raise UclustParseError on duplicate seed ID declarations
        """
        with open(self.uc_incorrect_1_fp, 'U') as uc_f:
            self.assertRaises(UclustParseError,
                              clusters_from_uc_file,
                              uc_f)


# VSEARCH test uc output
uc_output = """S    0   100 *   *   *   *   *   s1_80   *
H   0   100 100.0   *   0   0   *   s1_81   s1_80
H   0   100 100.0   *   0   0   *   s1_82   s1_80
S   1   100 *   *   *   *   *   s1_0    *
H   1   100 100.0   *   0   0   *   s1_1    s1_0
S   2   100 *   *   *   *   *   s1_10   *
H   2   100 100.0   *   0   0   *   s1_12   s1_10
H   2   100 100.0   *   0   0   *   s1_13   s1_10
S   3  100 *   *   *   *   *   s1_118  *
H   3  100 100.0   *   0   0   *   s1_119  s1_118
S   4  100 *   *   *   *   *   s1_11   *
S   5  100 *   *   *   *   *   s1_25   *
"""

# duplicate seed declarations (s1_80)
uc_output_incorrect_1 = """S    0   100 *   *   *   *   *   s1_80   *
H   0   100 100.0   *   0   0   *   s1_81   s1_80
H   0   100 100.0   *   0   0   *   s1_82   s1_80
S   1   100 *   *   *   *   *   s1_80    *
H   1   100 100.0   *   0   0   *   s1_1    s1_80
S   2   100 *   *   *   *   *   s1_10   *
H   2   100 100.0   *   0   0   *   s1_12   s1_10
H   2   100 100.0   *   0   0   *   s1_13   s1_10
S   3  100 *   *   *   *   *   s1_118  *
H   3  100 100.0   *   0   0   *   s1_119  s1_118
S   4  100 *   *   *   *   *   s1_11   *
S   5  100 *   *   *   *   *   s1_25   *
"""

# Test dereplicating sequences using default parameters
seqs_to_derep = """>HWI-ST157_0368:1:2102:15078:69955#0/1
TACGTAGGGCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGTGCGCAGGCGGTCTGTTAAGTCTGTAGTTAAAGGCTGTGGCTCAACTATGGTTAGTT
>HWI-ST157_0368:1:2103:7895:197066#0/1
TACGTAGGGGGCAAGCGTTGTCCGAATTTACTGGGTGTAAAGGGAGCGCAGACGGCACGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:1207:16180:126921#0/1
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGATACTTAAGTCTGGTGTGAAAACCTAGGGCTCAACCCTGGGACTGC
>HWI-ST157_0368:1:1106:11378:83198#0/1
TACGTAGGGAGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTTGCAAGTCAGATGTGAAATCTATGGGCTCAACCCATAAACTGC
>HWI-ST157_0368:1:2103:7895:197066#0/2
TACGTAGGGGGCAAGCGTTGTCCGAATTTACTGGGTGTAAAGGGAGCGCAGACGGCACGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:2103:7895:197066#0/3
TACGTAGGGGGCAAGCGTTGTCCGAATTTACTGGGTGTAAAGGGAGCGCAGACGGCACGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:1207:16180:126921#0/2
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGATACTTAAGTCTGGTGTGAAAACCTAGGGCTCAACCCTGGGACTGC
>HWI-ST157_0368:1:1207:16180:126921#0/3
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGATACTTAAGTCTGGTGTGAAAACCTAGGGCTCAACCCTGGGACTGC
"""

# Test dereplicating a file which is a concatenation of two separately
# dereplicated files. The input fasta file contains abundance information.
seqs_to_derep_merged_derep_files = """>HWI-ST157_0368:1:2102:15078:69955#0/1;size=1;
TACGTAGGGCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGTGCGCAGGCGGTCTGTTAAGTCTGTAGTTAAAGGCTGTGGCTCAACTATGGTTAGTT
>HWI-ST157_0368:1:2103:7895:197066#0/1;size=3;
TACGTAGGGGGCAAGCGTTGTCCGAATTTACTGGGTGTAAAGGGAGCGCAGACGGCACGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:1207:16180:126921#0/1;size=3;
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGATACTTAAGTCTGGTGTGAAAACCTAGGGCTCAACCCTGGGACTGC
>HWI-ST157_0368:1:1106:11378:83198#0/1;size=1;
TACGTAGGGAGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTTGCAAGTCAGATGTGAAATCTATGGGCTCAACCCATAAACTGC
>HWI-ST157_0368:1:2102:15078:69955#1/1;size=1;
TACGTAGGGCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGTGCGCAGGCGGTCTGTTAAGTCTGTAGTTAAAGGCTGTGGCTCAACTATGGTTAGTT
>HWI-ST157_0368:1:2103:7895:197066#1/1;size=3;
TACGTAGGGGGCAAGCGTTGTCCGAATTTACTGGGTGTAAAGGGAGCGCAGACGGCACGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:1207:16180:126921#1/1;size=3;
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGATACTTAAGTCTGGTGTGAAAACCTAGGGCTCAACCCTGGGACTGC
>HWI-ST157_0368:1:1106:11378:83198#1/1;size=1;
TACGTAGGGAGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTTGCAAGTCAGATGTGAAATCTATGGGCTCAACCCATAAACTGC
"""

# Sequences to dereplicate with final clusters as follows:
# 2 clusters with abundance 6 and 7
# 3 clusters with abundance 1
# 3 clusters with abundance 2, 3, 5
seqs_to_derep_max_min_abundance = """>HWI-ST157_0368:1:1106:10560:153880#0/1
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:1106:10560:153880#0/2
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:1106:10560:153880#0/3
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:1106:10560:153880#0/4
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:1106:10560:153880#0/5
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:1106:10560:153880#0/6
TACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGCGTAGGCGGGAATGCAAGTCAGATGTGAAATCCAGGGGCTTAACCCTTGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/1
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/2
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/3
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/4
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/5
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/6
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:2104:14337:180515#0/7
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGGGCGCGCAGGCGGTCTGGCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC
>HWI-ST157_0368:1:1102:8490:14349#0/1
AACGTAGGTCACAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGCAGGCGGGAAGACAAGTTGGAAGTGAAATCTATGGGCTCAACCCATAAACTGC
>HWI-ST157_0368:1:1205:18016:113727#0/1
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTGATTAAGTTAGATGTGAAATCCCCGGGCTTAACCTGGGGATGGC
>HWI-ST157_0368:1:1201:16382:127646#0/1
TACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCGGTTTCGTAAGTCGGGTGTGAAAGGCGGGGGCTTAACGCCCGGACTGG
>HWI-ST157_0368:1:2103:12440:90119#0/1
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCACGCCAAGTCAGCGGTGAAATTTCCGGGCTCAACCCGGAGTGTGC
>HWI-ST157_0368:1:2103:12440:90119#0/2
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCACGCCAAGTCAGCGGTGAAATTTCCGGGCTCAACCCGGAGTGTGC
>HWI-ST157_0368:1:1106:15269:103850#0/1
TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGGTTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGC
>HWI-ST157_0368:1:1106:15269:103850#0/2
TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGGTTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGC
>HWI-ST157_0368:1:1106:15269:103850#0/3
TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGGTTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGC
>HWI-ST157_0368:1:1205:9745:86166#0/1
TACGTAGGTCCCGAGCGTTATCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCCGGGCTCAACTTGGGTGCTGC
>HWI-ST157_0368:1:1205:9745:86166#0/2
TACGTAGGTCCCGAGCGTTATCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCCGGGCTCAACTTGGGTGCTGC
>HWI-ST157_0368:1:1205:9745:86166#0/3
TACGTAGGTCCCGAGCGTTATCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCCGGGCTCAACTTGGGTGCTGC
>HWI-ST157_0368:1:1205:9745:86166#0/4
TACGTAGGTCCCGAGCGTTATCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCCGGGCTCAACTTGGGTGCTGC
>HWI-ST157_0368:1:1205:9745:86166#0/5
TACGTAGGTCCCGAGCGTTATCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGATGATTAAGTGGGATGTGAAATACCCGGGCTCAACTTGGGTGCTGC
"""

# Test sort by abundance functionality in VSEARCH
seqs_to_sort = """>HWI-ST157_0368:1:2105:3428:36721#0/1;size=5;
TACGTAGGGTGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGCTCGTAGGCGGTTCGTCGCGTCCGGTGTGAAAGTCCATCGCTTAACGGTGGATCTGC
>HWI-ST157_0368:1:2106:18272:88408#0/1;size=2;
TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCGGAGCAAGTCTGAAGTGAAAGCCCGGGGCTCAACCCCGGGACTGC
>HWI-ST157_0368:1:1106:12200:200135#0/1;size=1;
TACGTAGGTGGCGAGCGTTATCCGGAATTATTGGGCGTAAAGAGGGAGCAGGCGGCACTAAGGGTCTGTGGTGAAAGATCGAAGCTTAACTTCGGTAAGCC
>HWI-ST157_0368:1:1201:8401:113582#0/1;size=10;
TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTGTGGCAAGTCTGATGTGAAAGGCATGGGCTTAACCTGTGGACTGC
>HWI-ST157_0368:1:2208:9135:145970#0/1;size=1;
TACGTAGGGGGCGAGCGTTGTCCGGAATTACTGGGCGTAAGGGGAGCGTAGGCGGTCGATTAAGTTAGATGTGAAACCCCCGGGCTTAACTTGGGGACTGC
>HWI-ST157_0368:1:2204:20491:181552#0/1;size=10;
TACGTAGGGTGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGCTCGTAGGCGGTTCGTCGAGTCTGGTGTGAAAGTCCATCGCTTAACGGTGGATCCGC
>HWI-ST157_0368:1:2105:6731:137157#0/1;size=4;
TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTGATAAGTCTGAAGTTAAAGGCTGTGGCTCAACCATAGTTCGCT
>HWI-ST157_0368:1:2107:19923:3944#0/1;size=100;
TGCATTTTCTCTTATCGAAAACCTTCAGCGTTCTGATCTGAATCCCGTCGAAGAGGCTAAGGGCTATCGCCAACTCATTGATGCCAGCGGGATGACCCAGG
"""

# Grinder simulated chimeric reads using Greengenes 13.8 release
# command: grinder -random_seed 100 -reference_file 97_otus_gg_13_8.fasta \
#                -forward_reverse ./primers.fna -length_bias 0 -copy_bias 1 \
#                -unidirectional 1 -read_dist 150 -mutation_dist uniform 0.1 \
#                -mutation_ratio 100 0 -total_reads 1000 -diversity 10 \
#                -chimera_perc 10 -od grinder_chimeric_reads_illumina
# primers.fna contain 515f and 806r primers
# reads >251 reference=4370324,646991 amplicon=488..779,499..789
#    >320 reference=646991,4370324,646991 amplicon=499..789,488..779,499..789
#    >36 reference=4370324,814974 amplicon=488..779,479..769
#    >672 reference=814974,160832 amplicon=479..769,436..727
#    >142 reference=160832,814974 amplicon=436..727,479..769 errors=2%G
#    >201 reference=4304512,510574 amplicon=451..742,501..793
#    >241 reference=646991,4370324 amplicon=499..789,488..779 errors=13%A
#    >279 reference=311922,160832,510574 amplicon=481..773,436..727,501..793
#    >299 reference=4370324,4304512 amplicon=488..779,451..742
#    >359 reference=646991,4370324 amplicon=499..789,488..779 errors=52%A
#    >375 reference=4304512,769294 amplicon=451..742,504..795
#    >407 reference=4304512,579954 amplicon=451..742,488..779
#    >423 reference=4370324,579954 amplicon=488..779,488..779
#    >516 reference=814974,579954 amplicon=479..769,488..779
#    >618 reference=814974,646991 amplicon=479..769,499..789 errors=32%C
#    >717 reference=814974,510574 amplicon=479..769,501..793
#    >902 reference=510574,579954 amplicon=501..793,488..779
#    >918 reference=814974,4370324 amplicon=479..769,488..779
#    >941 reference=579954,4304512 amplicon=488..779,451..742
# are chimeric
amplicon_reads = """>3;size=102;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>16;size=95;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>22;size=93;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>2;size=87;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>39;size=84;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>4;size=79;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>6;size=72;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>11;size=70;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>45;size=67;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>1;size=65;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>251;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>30;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>320;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>36;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>425;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>672;size=2;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>10;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GGGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>100;size=1;
TTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>102;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGACTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>115;size=1;
GTGCCAGCAGCCGCGGTATTACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>123;size=1;
GTGACAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>132;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGACACTGCAAGTCTTGAGATCGGAAG
>134;size=1;
GTGCCGGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>14;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCACGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>140;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGGAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>142;size=1;
GGGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>144;size=1;
GTGTCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>148;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTATCATTCTTGAGTATAGATG
>15;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGAGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>156;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGACGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>161;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGGGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>162;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGCTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>186;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTAATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>201;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>203;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAACGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>21;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTCAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>217;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGTAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>218;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAACGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>221;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAATACTGGAG
>222;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGGGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>225;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGAGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>233;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTGTACTTGAGTGTTGTAA
>234;size=1;
GTGCCAGCAGCCTCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>235;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAACGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>24;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAAGTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>241;size=1;
GTGCCAGCAGCCACGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGTCATTCTTGAGTATAGATG
>249;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTCAAACTGCAAGTCTTGAGATCGGAAG
>259;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCGCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>26;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACATGATACTGCCTTGCTCGAGTACTGGAG
>266;size=1;
GTGCCAGCGGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>27;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGTGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>279;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>296;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGCGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>299;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>303;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCAGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>306;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAA
>307;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGTTACTGGTATACTTGAGTGTTGTAA
>32;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGGTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>322;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGGCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>326;size=1;
GTGCCAGCAGCCGCGGTAATACGGTGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>332;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGGTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>333;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTATGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>338;size=1;
GTGCCAGCAGCCGCGGTATTACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>359;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTAGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGTCATTCTTGAGTATAGATG
>360;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGTTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>362;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGAATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>364;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCATTTAGAACTGGTTAACTAGAGTATTGGAG
>366;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCCTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>369;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTAATGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>37;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACTTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>371;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGCGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>373;size=1;
GTGCCAGCAGACGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>374;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGCTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>375;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>386;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGTAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>387;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACAGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>392;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAACCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>393;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGTGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>397;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTTCCATTGATACTGGTATACTTGAGTGTTGTAA
>405;size=1;
GTACCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>407;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>414;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGTGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>418;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTTAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>423;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>431;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AGGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>436;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGTAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>444;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTGAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>445;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGCGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>456;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGCCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>460;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
CAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>469;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGAGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>470;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCGCGAGTACTGGAG
>477;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCGGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>479;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGGGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>486;size=1;
GTGCCAGCAGCTGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>500;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTACAAGTCTTGAGATCGGAAG
>515;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGCAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>516;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>528;size=1;
GTTCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>530;size=1;
GTGCCAGCAGGCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>531;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGGCTTGTAG
>549;size=1;
GTCCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>551;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCAGTTGATACTGCCTTGCTCGAGTACTGGAG
>557;size=1;
GTGCCAGCAGCCGCTGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>559;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGCTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>561;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCATGTCTTGAGATCGGAAG
>562;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCATTGATACTGGATGTCTTGAGTGTGAGAG
>564;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAATTGCCATTGATACTGTCATTCTTGAGTATAGATG
>566;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTAATGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>568;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTTAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>57;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGAGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>570;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGCGTTACATAGA
>578;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATGTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>586;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATATCTTGAGTGTGAGAG
>587;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>596;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGTGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>600;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGTTTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>612;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCAGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>618;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACCAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>625;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTATTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>632;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGCGCAGACTTGAGTGATGTAG
>649;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCGACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>650;size=1;
GTGCCAGCAGCCGTGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>651;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACCGGCAGACTTGAGTGATGTAG
>66;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTGGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>664;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTGATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>673;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGCGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>675;size=1;
GTGCCAGCAGCCGCGGTAATACCTAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>682;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTTT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>690;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTGGATACTGCCTTGCTCGAGTACTGGAG
>699;size=1;
GTTCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>7;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCATAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>709;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTCTGGCTTGAGTTCGGCAG
>717;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>73;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCGCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>740;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGTGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>745;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGAGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>746;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGTCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>748;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCATTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>760;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGCGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>766;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTAATACTGGTATACTTGAGTGTTGTAA
>77;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCGGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>778;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTTGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>791;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
TAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>797;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTTTGGCTTGAGTTCGGCAG
>809;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGAGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>813;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACAGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>814;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>816;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGAGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>817;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>82;size=1;
GTGCCAGCAGGCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>821;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGTGTACTGGAG
>824;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGTGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>827;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTAGTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>83;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAAGTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>842;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGGGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>851;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGTCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>853;size=1;
GTTCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>862;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTGGAGTGTTGTAA
>863;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTGAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>866;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAG
>871;size=1;
GCGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>879;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCTGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>886;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGAGG
>892;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTT
AAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGGAGGCTAGAGTCTTGTAG
>895;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGATCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>897;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGAATTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>902;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>904;size=1;
CTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>91;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAC
>912;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCATTTTAGAACTGGTTAACTAGAGTATTGGAG
>916;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGGGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>918;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>920;size=1;
GTGCCAGCAGCCGCGGTAATACGTCGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>921;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGCGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>925;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGGAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>930;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGT
AAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAG
>941;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAG
>942;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTATAGTATTGGAG
>945;size=1;
GTGCCAGCAGCCGTGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
>947;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTCATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>948;size=1;
GTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTT
AAGTCTAGAGTTAAAACCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGA
>952;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGCGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>956;size=1;
GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTT
AAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCGTTTGATACTGGATGTCTTGAGTGTGAGAG
>958;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCCGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTC
GCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAG
>964;size=1;
GTGCCAGCAGCCGCCGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATT
AAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAG
>967;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTAGAG
>984;size=1;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTAGAGATG
>992;size=1;
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGT
AAATCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAG
>993;size=1;
GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGACTACAT
AAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAG
"""

# Single chimeric sequence (251) with both parents to test alignment
# and tabular output de novo
single_chimera = """>22;size=93;
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATT
AAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAA
>45;size=67;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATG
>251;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGGTATACTTGAGTGTTGTAA
"""

# Alignment for chimeric sequence (251) against parents (22 and 45) de novo
single_chimera_aln = """
------------------------------------------------------------------------
Query   (  150 nt) 251;size=2;
ParentA (  150 nt) 45;size=67;
ParentB (  150 nt) 22;size=93;

A     1 GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT 80
Q     1 GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT 80
B     1 GTGCCAGCAGCCGCGGTAATACGGAGGgTGCGAGCGTTgTCCGGATTTATTGGGTTTAAAGGGTaCGTAGGCGGtgTatT 80
Diffs                              A          A                         A         AA AA 
Votes                              +          +                         +         ++ ++ 
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A    81 AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGtcATtCTTGAGTaTaGatg 150
Q    81 AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGGTATACTTGAGTGTTGTAA 150
B    81 AAGTCAGTGGTGAAAgCCTGCgGCTCAACcGTAGgagTGCCATTGATACTGGTATACTTGAGTGTTGTAA 150
Diffs                  A     A       A    AAA              BB  B       B B BBB
Votes                  +     +       +    +++              ++  +       + + +++
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAxxxxxxxxxxxxxxBBBBBBBBBBBBBBBBBBB

Ids.  QA 94.7%, QB 91.3%, AB 86.0%, QModel 100.0%, Div. +5.6%
Diffs Left 13: N 0, A 0, Y 13 (100.0%); Right 8: N 0, A 0, Y 8 (100.0%), Score 0.8291
"""

# Tabular format for UCHIME output
single_chimera_tab = """0.0000\t22;size=93;\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\t45;size=67;\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.8291\t251;size=2;\t45;size=67;\t22;size=93;\t45;size=67;\t100.0\t94.7\t91.3\t86.0\t94.7\t13\t0\t0\t8\t0\t0\t5.3\tY
"""

# Single chimeric sequence for reference chimera checking
single_chimera_ref = """>251;size=2;
GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGT
AAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGGTATACTTGAGTGTTGTAA
"""

# Reference database for UCHIME ref
uchime_ref_db = """>4304512
TGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGCGTCCTTCGGGACGAGTGGCAGACGGGTGAGTAACGCGTGGGAACGTACCCTTTGGTTCGGAACAACTCCGGGAAACTGGAGCTAATACCGGATAAGCCCTTCGGGGGAAAGATTTATCGCCTTTAGAGCGGCCCGCGTCTGATTAGCTAGTTGGTGGTGTAATGGACCACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGACCAGCCACATTGGGACTGAGACACGGCTCAAACTCCTACGGGAGGCAGCAGTGGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCCATGCCGCGTGTATGATGAAGGTCTTAGGATTGTAAAATACTTTCACCGGTGAAGATAATGACTGTAGCCGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTTAAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGATGTCTTGAGTGTGAGAGAGGTATGTGGAACTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACATACTGGCTCATTACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATTGCTAGTTGTCGGGATGCATGCATTTCGGTGACGCAGCTAACGCATTAAGCAATCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCACCTTTTGACATGCCTGGACCGCCAGAGAGATCTGGCTTTCCCTTCGGGGACTAGGACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCATTAGTTGCCATCATTTAGTTGGGAACTCTAATGGGACTGCCGGTGCTAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACAGGGTGGGCTACACACGTGCTACAATGGCGACTACAGAGGGTTAATCCTTAAAAGTCGTCTCAGTTCGGATTGTCCTCTGCAACTCGAGGGCATGAAGTTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGTTCTACCCGAAGGCGCTGCGCTAACCGCAAGGGGGCAGGCGACCACGGTAGGGTCAGCGACTGGGGTGAAGTCGTAACAAGGTAACC
>814974
GATGAACGCTAGCCGTGTGCCTAATACATGCATGTCGTACGAGAGTACTTGTACTCTAGTGGCGAATGGGTGAGTAACACGTACCTAACCTACTTTTAAGATTGGAATAACTACTGGAAACAGTAGCTAATGCCGAATACGTATTAACTTCGCATGAAGATAATATAAAAGGAGCGTTTGCTCCGCTTAGAAATGGGGGTGCGCCATATTAGTTAGTTGGTAGGGTAATGGCCTACCAAGACGATGATATGTAGCCGGGCTGAGAAGCTGATCGGCCACACTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATTTTCCGCAATGAGCGAAAGCTTGACGGAGCGACACGGCGTGCAGGATGAAGGTCTTCGGATCGTAAACTGCTGTGGTTAGGGAAGAAAAGCAAAATAGGAAATGATTTTGCCCTGACGGTACCTAACTAGAAAGTGACGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGTCACAAGCGTTATCCGGATTTATTGGGCGTAAAGCGTTCGTAGGCGGTTTGTTAAGTCTAGAGTTAAAGCCTGGGGTTCAACCCCAGCCCGCTTTGGATACTGACAAACTAGAGTTACATAGAGGTAAGCGGAATTCTTAGTGAAGCGGTGGAATGCGTAGATATTAGGAAGAACACCAATGGCGAAGGCAGCTTACTGGATGTACACTGACGCTCAGGGACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATCACTAGCCGCTAGAAAATTTAGTGGCACAGCTAACGCATTAAGTGATCCGCCTGAGTAGTATGCTCGCAAGAGTGAAACTTAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTTGAAGATACGCGTAGAACCTTACCCACTCTTGACATCTTCCGCAATGCTACAGAGATGTAGTGGAGGTTAACGGAATGACAGATGGTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTTTGGTTAAGTCCAGCAACGAGCGCAACCCTTGTCTTTAGTTACTAATATTAAGTTAAGGACTCTAGAGAGACTGCCAGGGTAACCTGGAGGAAGGTGGGGACGACGTCAAATCATCATGCCTCTTACGAGTGGGGCAACACACGTGCTACAATGGTCGGTACAAAGAGAAGCAAGATGGCGACATGGAGCAAACCTCAAAAAACCGATCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGAAGTTGGAATCGCTAGTAATCGTAGATCAGCTACGCTACGGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCACACCATGGGAGCTGGTAATGCCCGAAGCCGGTTAGTTAACTTCGGAGACGACTGTCTAAGGCAGGACTGGTGACTGGGGTG
>160832
GTCGAGCGGCGGACGGGTGAGTAACGGCTGGGAACCTGCCCTGACGCGGGGGATAACCGTTGGAAACGACGGCTAATACCGCATAATGTCTTAGTTCATTACGAGCTGGGACCAAAGGTGGCCTCTACATGTAAGCTATCGCGTTGGGATGGGCCCAGTTAGGATTAGCTAGTTGGTAAGGTAATGGCTTACCAAGGCRACGATCCTTAKCTGGTTTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGGGAGACCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGCAGTGAGGAAGGTGGTGTACTTAATAAGTGCATGGCTTGACGTTAGCTGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGCGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAACCGCATTTGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCAGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGATGCGAAAGCGTGGGTAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGCTGTCTACTTGGAGGTTGAGGTTTAAGACTTTGGCTTTCGGCGCTAACGCATTAAGTAAACCGCCGGGGGAGTACGGCCGCAGGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGGGGTTTAATTCGATGCAACGCGAAAAACCTTACCTACTCTTGACATCCAACGAATCCTTTAAAGATGGATGAGTGCCTTCGGAAGCGCTGAAACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTGTTTGCCAGCACATAATGGTGGGAACTCCAGGGAGACTGCCGGTGATAAACCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGAGTAGGGCTACACACGTGCTACAATGGCAGATACAGAGGGCAGCGAAGCTGCGAAGTGGAGCGAATCCCTTAAAGTCTGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGCTGCACCAGAAGTAGATAGCTTAACCTTCGGGAGGGCGTTTACCACGGTGTGGTTCATGACTGGGGTGAAGTCGTAACAAGGTAGCCCTAGGGGAACCTGCGGCTG
>573757
GCCTAAGGCAGGCAAGTCGAACGATGATCTCCAGCCTGCTGGGGGGATTAGAGGAGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTCTGGGATAAGCCTGGGAAACTGGGTCTAATACTGGATACGACCTTCCCACGCATGTGGTGTTGGTGGAAAGCTTTTGTGGTTTTGGATGGACTCGCGGCCTATCAGCTTGTTGGTGGGGTAATGGCCTACCAAGGCGACGACGGGTAGCCGGCCTGAGAGGGTGGACGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGAGGGATGAAGGCCTTCGGGTTGTAAACCTCTTTCAGTAGGGAAGAAGCGAAAGTGACGGTACCTGCAGAAGAAGCGCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTCGCGTCTGCCGTGAAAGTCCGGGGCTTAACTCCGGATCTGCGGTGGGTACGGGCAGACTTGAGTGATGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCATTAACTGACGCTGAGGAGCGAAAGCATGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACACCGTAAACGTTGGGCGCTAGGTGTGGGGCCTATTCCATGGGTTCCGTGCCGCAGCTAACGCATTAAGCGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGTTTGACATATGCCGGAAAGCCGTAGAGATACGGCCCCTTTTTGTCGGTATACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTCCTATGTTGCCAGCACGCCCGTAGGGTGGTGGGGACTCATAGGAGGCTGCCGGGGTCAACTCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTATGTCCAGGGCTTCACGCATGCTACAATGGCCGGTACAAAGGGCTGCGATCCCGTGAGGGGGAGCGAATCCCAAAAAGCCGGTCTCAGTTCGGATTGGGGTCTGCAACTCGACCCCATGAAGTCGGAGTCGCTAGTAATCGCAGATCAGCAACGCTGCGGTGAATACGTTCCCGGGCCTTGCACACACCGCCCGTCAC
>579954
TGGCGGCGTGGATAAGACATGCAAGTCGAACGGGATATTGTTTGTAGCAATACAAGCGATGTCTAGTGGCGTAAGGGTGCGTAACACGTGGGGAATCTGCCGAGAAGTGGGGGATAGCTCGCCGAAAGGCGAATTAATACCGCATGTGGTTAGGGAAGACATCTTCCCGACACTAAAGCCGGGGCAACCTGGCGCTTCTTGATGACCCCGCGGCCTATCAGCTAGTCGGTGAGGTAACGGCTCACCAAGGCTATGACGGGTAGCTGGTCTGAGAGGACGACCAGCCACACTGGAACTGAGACACGGTCCAGACACCTACGGGTGGCAGCAGTCGAGAATTTTTCTCAATGGGGGAAACCCTGAAGGAGCGACGCCGCGTGGAGGATGAAGGTCTTCGGATTGTAAACTCCTGTCATGCGGGAACAATTGTCACCGATTAACTGTCGGGGGCTTGATAGTACCAGAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCGGAGCTCAACTCCGAAACTGCACTTGATACTGCCTTGCTCGAGTACTGGAGAGGAGATTGGAATTTACGGTGTAGCAGTGAAATGCGTAGATATCGTAAGGAAGACCAGTGGCGAAGGCGAATCTCTGGACAGTTACTGACACTGAGGCACGAAGGCCAGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCTGGCAGTAAACGGTGCGCGTTTGGTGTGGGAGGATTCGACCCCTTCTGTGCCGGAGCTAACGCGTTAAACGCGCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGAAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGCTTAATTCGATGCAACGCGAAGAACCTTACCTAGCCTTGACATGCATCTCTAAGTCGGTGAAAGCCGGCGACTATAGCAATATAGAATTTGCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGTGAACTGTTGCCACCGATCTTCGGATCGAGCACTCTGTTCAGACTGCCCTGTGAAACGGGGAGGAAGGTGGGGATGACGTCAAGTCAGCATGGCCCTTACGGCTAGGGCTGCACACGTACTACAATGCTCAGTACAGAATGAACCGAAACCGCGAGGTGGAGGAAATCCATAAAACTGAGCCCAATTCGGATTGGAGGCTGCAACTCGCCTCCATGAAGCCGGAATCGCTAGTAATGGCGTATCAGCTACGACGCCGTGAATACGTTCCCGGGCCTTGTACACACCG
>311922
GATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACATGAAGTGCTTGCACTTTGATGACGAGTGGCGGACGGGTGAGTAATGCTTGGGAATTTGCCTTTGCGCGGGGGATAACCATTGGAAACGATGGCTAATACCGCATAATGTCTACGGACCAAAGGGGGCTTAGGCTCCCACGTGAAGAGAAGCCCAAGTGAGATTAGCTAGTTGGTGGGGTAAAGGCTCACCAAGGCGACGATCTCTAGCTGTTCCGAGAGGAAGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCGCAATGGGGGGAACCCTGACGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGTTATGAGGAAAGGTTGTTGGTTAATACCCAGCAGCTGTGACGTTAATAACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTAATTAAGCTAGAAGTGAAAGCCCTGCGCTCAACGTGGGAAGGCCTTTTAGAACTGGTTAACTAGAGTATTGGAGAGGGGAGTGGAATTCCAGGTGTAGCGGTGAAATGCGTAGATATCTGGAGGAACACCGGTGGCGAAGGCGACTCTCTGGCCCAAATACTGACGCTCATGTGCGAAAGTGTGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCACACCGTAAACGATGTCTACTAGCTGTGTGCGAATGTATATTTGTGCGTAGCGCAGCCAACGCGATAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTACACTTGACATGCAGAGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACTCTGACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCTTATTTGCCAGCATTAAGTTGGGGACTTTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGTGTAGGGCTACACACGTGCTACAATGGTAAGTACAGAGGGAAGCGAACTTGTGAGAGTAAGCGGACCCCAAAAAGCTTATCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGATGCAAAAGAAGTAGGTAGCATAACCTTTAGGAGTGCGCTTACCACTTTGTGTTTCATGACTGGGGT
>4370324
TCAGGATGAACGCTAGCGACAGGCCTAACACATGCAAGTCGAGGGGTAACATTGGTAGCTTGCTACCAGATGACGACCGGCGCACGGGTGAGTAACGCGTATGCAACCTTCCTTTAACAGGAGAATAGCCCCCGGAAACGGGGATTAATGCTCCATGGCACTCTAATTTCGCATGGAATAAGAGTTAAAGTTCCGACGGTTAAAGATGGGCATGCGTGACATTAGCCAGTTGGCGGGGTAACGGCCCACCAAAGCAACGATGTCTAGGGGTTCTGAGAGGAAGGTCCCCCACACTGGTACTGAGACACGGACCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGTCGCGTGCAGGATGACTGCCCTATGGGTTGTAAACTGCTTTTGTACGGGAAGAAATGTACTTACGAGTAAGTATTTGCCGGTACCGTACGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGTAAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATGAGGTAGGCGGAATGAGTAGTGTAGCGGTGAAATGCATAGATATTACTCAGAACACCAATTGCGAAGGCAGCTTACTAAACTATAACTGACGCTGAAGCACGAAAGCGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATTACTGGTTGTTTGCAATACACCGCAAGCGACTGAGCGAAAGCATTAAGTAATCCACCTGGGGAGTACGTCGGCAACGATGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGTTTAAATGGGAAGTGACAGGGGTAGAAATACCTTTTTCTTCGGACACTTTTCAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGGTTAAGTCCCATAACGAGCGCAACCCCTGTTGTTAGTTACCAGCATGTAAAGATGGGGACTCTAACAAGACTGCCGGTGTAAACCGCGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACATCCAGGGCTACACACGTGTTACAATGGCAGGTACAAAGGGCAGCTACACAGCGATGTGATGCTAATCTCGAAAACCTGTCCCAGTTCGGATTGAAGTCTGCAACCCGACTTCATGAAGCTGGAATCGCTAGTAATCGCGCATCAGCCATGGCGCGGTGAATACGTTCCCGGGCCTTGTACACTCCGCCCGTCAAGCCATGGAAGCCGGGAGTACCTGAAG
>646991
AGTTTGATCTTGGCTCAGGATGAACGCTAGCGGCAGGCCTAATACATGCAAGTCGTGGGGCATCAGCGCCTTCGGGCGGCTGGCGACCGGCGCACGGGTGCGTAACGCGTATGCAACCTGCCCACAACAGGGGGACAGCCTTCGGAAACGAGGATTAATACCCCATGATACAGGGGTACCGCATGGTGCCTTTCGTCAAAGGTTTCGGCCGGTTGTGGATGGGCATGCGTCCCATTAGCTAGTAGGCGGGGTAACGGCCCACCTAGGCTATGATGGGTAGGGGTTCTGAGAGGACGATCCCCCACACTGGTACTGAGATACGGACCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGGGCAATGGGCGGAAGCCTGACCCAGCCATGCCGCGTGCAGGACGAAGGCCCTCGGGTCGTAAACTGCTTTTATACGGGAAGAACTGCGTCCTGCGGGACGCGCTGACGGTACCGTACGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATTAAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAAGGGTGGGCGGAATTCCGCATGTAGCGGTGAAATGCATAGATATGCGGAGGAACACCGAGAGCGAAGGCAGCTCACTAGGCACGACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGTAACTAGGTGTGTGCGACACAGAGTGCGCGCCCAAGCGAAAGCGATAAGTTACCCACCTGGGGAGTACGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGGGTCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGCTCGAATGGCCTATGACAGGCCCAGAGATGGGCCCTTCCTCGGACATAGGTCAAGGTGCTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCCCTAGTTGCCATCAGGTAAAGCTGGGGACTCTAGGGGGACTGCCTGCGCAAGCAGAGAGGAAGGAGGGGACGATGTCAAGTCATCATGGCCCTTACGCCCAGGGCTACACACGTGCTACAATGGCGCATACAGAGGGTAGCCACCTGGCGACAGGGCGCCAATCTCAAAAAGTGCGTCTCAGTTCGGATCGGGGCCTGCAACTCGGCCCCGTGAAGTCGGAATCGCTAGTAATCGCAGATCAGCCATGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGAAGCCGGGGGCACCTGAAGTCGGGGGTAACAACCCGCCTAGGGTGAAACTGGTAACTGGGGCTAAGTCGTAACAAGGTAACCGTA
>510574
GAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGAATGCTTTACACATGCAAGTCGAGCGGCAGCGCGGGGGCAACCCTGGCGGCGAGCGGCGAACGGGTGAGTAACACATCGGAACGTACCCAATTGAGGGGGATAGCCCGGCGAAAGCCGGATTAATACCGCATAAGTCCTGAGGGAGAAAGCGGGGGACCGCAAGGCCTCGCGCGATTGGAGCGGCCGATGTCGGATTAGCTAGTTGGTGGGGTAAAGGCTCACCAAGGCGACGATCCGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCCATTCCGCGTGAGTGAAGAAGGCCTTCGGGTTGTAAAGCTCTTTCGGACGGAAAGAAATCGCCCGGGTAAATAATCCGGGTGGATGACGGTACCGTAAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTACATAAGACAGGTGTGAAATCCCCGGGCTCAACCTGGGAATGGCGCTTGTGACTGTGTGGCTTGAGTTCGGCAGAGGGGGGTGGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAACACCGATGGCGAAGGCAGCCCCCCTGGGCCGTGACTGACGCTTATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCGACTAGGTGTTGGGGAAGGAGACTTCTTTAGTACCGTAGCTAACGCGTGAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAACGCGAAAAACCTTACCTACCCTTGACATGCCAGGAACCTTGCTGAGAGGTGAGGGTGCCCGAAAGGGAACCTGGACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAATTGCCATCATTGAGTTGGGCACTTTAATGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGGGTAGGGCTTCACACGTCATACAATGGCGCATACAGAGGGTTGCCAACCCGCGAGGGGGAGCGAATCCCAGAAAATGCGTCGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGGATCAGCATGTCGCGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTCCACCAGAAGTAGGTAGCTTAACCGCAAGGAGGGCGCTTACCACGGTGAGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACC
>769294
AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCAAGGGGAAAGTTTTCTTCGGAGAATTAGTATACTGGCGCACGGGTGAGTAATGTATAAGTAATCTACCTATAGGAAAGGAATAACTCTAAGAAATTGGGGCTAATACCATATAATGCAGCGGCACCGCATGGTGATGTTGTTAAAGTAATTTATTACGCCTATAGATGAGCTTGTATTCGATTAGCTTGTTGGTAAGGTAACGGCTTACCAAGGCGACGATCGATAGCTGGTCTGAGAGGATGATCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGGCAATGGACGAAAGTCTGACCCAGCAACGCCGCGTGGAGGATGAAGGTCGTAAGATCGTAAACTCCTTTTTTGGGGGAAGAAAAAACAGGTTTGTAGCCTGTATTGACTGTACCCTAAGAATAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGATTTACTGGGTATAAAGGGCTCGCAGGCGGGTTTGTAAGTCAGAGGTGAAATCCTACAGCTTAACTGTAGAACTGCCTTTGAAACTGCAAGTCTTGAGATCGGAAGAGAGAGATGGAATTCCAGGTGTAGTAGTGAAATACGTAGATATCTGGAAGAACACCAACGGCGAAGGCAGTCTCTTGGTCCGTATCTGACGCTGAGGAGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATACTAGGTGTCGGATTCTATGAATTCGGTGCCGCAGCTAATGCATTAAGTATTCCACCTGGGGAGTACGATCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCAGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTAGGCTTGAAATGCTGTGGACCGCTTGTGAAAGCAAGCTTCTCTTCGGAGCCGCAGTACAGGTGCTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCATTAGTTGCCATCAGGTTAAGCTGGGCACTCTAATGGGACTGCCTACGCAAGTAGTGAGGAAGGTGGGGATGACGTCAAGTCAGCATGGCCCTTACGCCTAGGGCTACACACGTGCTACAATGGATACTACAATGGGTTGCCAAGCCGCGAGGTGGAGCCAATCCCTTAAAAGTATCCTCAGTTCGGATTGGAGTCTGCAACTCGACTCCATGAAGCCGGAATTGCTAGTAATCGCGTATCAGCATGACGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGAAGCCGGGGGTACCCAAAGTCAGCGACCCAACCGCAAGGAGGGAGCTGCCTAAGGTAAAACTAGTGACTGGGGCTAAGTCGTAACAAGGTAACC
"""

# Reference database for single chimera sequence (ref)
uchime_single_ref_db = """>646991
AGTTTGATCTTGGCTCAGGATGAACGCTAGCGGCAGGCCTAATACATGCAAGTCGTGGGGCATCAGCGCCTTCGGGCGGCTGGCGACCGGCGCACGGGTGCGTAACGCGTATGCAACCTGCCCACAACAGGGGGACAGCCTTCGGAAACGAGGATTAATACCCCATGATACAGGGGTACCGCATGGTGCCTTTCGTCAAAGGTTTCGGCCGGTTGTGGATGGGCATGCGTCCCATTAGCTAGTAGGCGGGGTAACGGCCCACCTAGGCTATGATGGGTAGGGGTTCTGAGAGGACGATCCCCCACACTGGTACTGAGATACGGACCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGGGCAATGGGCGGAAGCCTGACCCAGCCATGCCGCGTGCAGGACGAAGGCCCTCGGGTCGTAAACTGCTTTTATACGGGAAGAACTGCGTCCTGCGGGACGCGCTGACGGTACCGTACGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTACGTAGGCGGTGTATTAAGTCAGTGGTGAAAGCCTGCGGCTCAACCGTAGGAGTGCCATTGATACTGGTATACTTGAGTGTTGTAAGGGTGGGCGGAATTCCGCATGTAGCGGTGAAATGCATAGATATGCGGAGGAACACCGAGAGCGAAGGCAGCTCACTAGGCACGACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGTAACTAGGTGTGTGCGACACAGAGTGCGCGCCCAAGCGAAAGCGATAAGTTACCCACCTGGGGAGTACGCTCGCAAGAGTGAAACTCAAAGGAATTGACGGGGGTCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGCTCGAATGGCCTATGACAGGCCCAGAGATGGGCCCTTCCTCGGACATAGGTCAAGGTGCTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCCCTAGTTGCCATCAGGTAAAGCTGGGGACTCTAGGGGGACTGCCTGCGCAAGCAGAGAGGAAGGAGGGGACGATGTCAAGTCATCATGGCCCTTACGCCCAGGGCTACACACGTGCTACAATGGCGCATACAGAGGGTAGCCACCTGGCGACAGGGCGCCAATCTCAAAAAGTGCGTCTCAGTTCGGATCGGGGCCTGCAACTCGGCCCCGTGAAGTCGGAATCGCTAGTAATCGCAGATCAGCCATGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGAAGCCGGGGGCACCTGAAGTCGGGGGTAACAACCCGCCTAGGGTGAAACTGGTAACTGGGGCTAAGTCGTAACAAGGTAACCGTA
>4370324
TCAGGATGAACGCTAGCGACAGGCCTAACACATGCAAGTCGAGGGGTAACATTGGTAGCTTGCTACCAGATGACGACCGGCGCACGGGTGAGTAACGCGTATGCAACCTTCCTTTAACAGGAGAATAGCCCCCGGAAACGGGGATTAATGCTCCATGGCACTCTAATTTCGCATGGAATAAGAGTTAAAGTTCCGACGGTTAAAGATGGGCATGCGTGACATTAGCCAGTTGGCGGGGTAACGGCCCACCAAAGCAACGATGTCTAGGGGTTCTGAGAGGAAGGTCCCCCACACTGGTACTGAGACACGGACCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGTCGCGTGCAGGATGACTGCCCTATGGGTTGTAAACTGCTTTTGTACGGGAAGAAATGTACTTACGAGTAAGTATTTGCCGGTACCGTACGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAATGGTAAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGTCATTCTTGAGTATAGATGAGGTAGGCGGAATGAGTAGTGTAGCGGTGAAATGCATAGATATTACTCAGAACACCAATTGCGAAGGCAGCTTACTAAACTATAACTGACGCTGAAGCACGAAAGCGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATTACTGGTTGTTTGCAATACACCGCAAGCGACTGAGCGAAAGCATTAAGTAATCCACCTGGGGAGTACGTCGGCAACGATGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCTGGGTTTAAATGGGAAGTGACAGGGGTAGAAATACCTTTTTCTTCGGACACTTTTCAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGGTTAAGTCCCATAACGAGCGCAACCCCTGTTGTTAGTTACCAGCATGTAAAGATGGGGACTCTAACAAGACTGCCGGTGTAAACCGCGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACATCCAGGGCTACACACGTGTTACAATGGCAGGTACAAAGGGCAGCTACACAGCGATGTGATGCTAATCTCGAAAACCTGTCCCAGTTCGGATTGAAGTCTGCAACCCGACTTCATGAAGCTGGAATCGCTAGTAATCGCGCATCAGCCATGGCGCGGTGAATACGTTCCCGGGCCTTGTACACTCCGCCCGTCAAGCCATGGAAGCCGGGAGTACCTGAAG
"""

# 3-way alignment for single chimeric sequence against reference
# database using UCHIME
single_chimera_ref_aln = """
------------------------------------------------------------------------
Query   (  150 nt) 251;size=2;
ParentA ( 1403 nt) 4370324
ParentB ( 1480 nt) 646991

A     1 tcaggatgaacgctagcgacaggcctaacacatgcaagtcgaggggtaacattggtagcttgctaccagatgacgaccgg 80
Q     1 -------------------------------------------------------------------------------- 0
B     1 agtttgatcttggctcaggatgaacgctagcggcaggcctaatacatgcaagtcgtggggcatcagcgccttcgggcggc 80
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A    81 cgcacgggtgagtaacgcgtatgcaaccttcctttaacaggagaatagcccccggaaacggggattaatgctccatggca 160
Q     1 -------------------------------------------------------------------------------- 0
B    81 tggcgaccggcgcacgggtgcgtaacgcgtatgcaacctgcccacaacagggggacagccttcggaaacgaggattaata 160
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   161 ctctaatttcgcatggaataagagttaaagttccgacggttaaagatgggcatgcgtgacattagccagttggcggggta 240
Q     1 -------------------------------------------------------------------------------- 0
B   161 ccccatgatacaggggtaccgcatggtgcctttcgtcaaaggtttcggccggttgtggatgggcatgcgtcccattagct 240
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   241 acggcccaccaaagcaacgatgtctaggggttctgagaggaaggtcccccacactggtactgagacacggaccagactcc 320
Q     1 -------------------------------------------------------------------------------- 0
B   241 agtaggcggggtaacggcccacctaggctatgatgggtaggggttctgagaggacgatcccccacactggtactgagata 320
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   321 tacgggaggcagcagtgaggaatattggtcaatggacgcaagtctgaaccagccatgtcgcgtgcaggatgactgcccta 400
Q     1 -------------------------------------------------------------------------------- 0
B   321 cggaccagactcctacgggaggcagcagtagggaatattgggcaatgggcggaagcctgacccagccatgccgcgtgcag 400
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   401 tgggttgtaaactgcttttgtacgggaagaaatgtacttacgagtaagtatttgccggtaccgtacgaataagcatcggc 480
Q     1 -------------------------------------------------------------------------------- 0
B   401 gacgaaggccctcgggtcgtaaactgcttttatacgggaagaactgcgtcctgcgggacgcgctgacggtaccgtacgaa 480
Diffs                                                                                   
Votes                                                                                   
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   481 taactcc-----------GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGG 549
Q     1 ------------------GTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGG 62
B   481 taagcaccggctaactccGTGCCAGCAGCCGCGGTAATACGGAGGgTGCGAGCGTTgTCCGGATTTATTGGGTTTAAAGG 560
Diffs                                                A          A                       
Votes                                                +          +                       
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

A   550 GTGCGTAGGCGGAATGGTAAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGtcATtCTTGAG 629
Q    63 GTGCGTAGGCGGAATGGTAAGTCAGTGGTGAAATCCTGCAGCTCAACTGTAGAGTTGCCATTGATACTGGTATACTTGAG 142
B   561 GTaCGTAGGCGGtgTatTAAGTCAGTGGTGAAAgCCTGCgGCTCAACcGTAGgagTGCCATTGATACTGGTATACTTGAG 640
Diffs     A         AA AA                A     A       A    AAA              BB  B      
Votes     +         ++ ++                +     +       +    +++              ++  +      
Model   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAxxxxxxxxxxxxxxBBBBBBBBBBB

A   630 TaTaGatgaggtaggcggaatgagtagtgtagcggtgaaatgcatagatattactcagaacaccaattgcgaaggcagct 709
Q   143 TGTTGTAA------------------------------------------------------------------------ 150
B   641 TGTTGTAAgggtgggcggaattccgcatgtagcggtgaaatgcatagatatgcggaggaacaccgagagcgaaggcagct 720
Diffs    B B BBB                                                                        
Votes    + + ++                                                                         
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A   710 tactaaactataactgacgctgaagcacgaaagcgtgggtatcaaacaggattagataccctggtagtccacgccgtaaa 789
Q   151 -------------------------------------------------------------------------------- 150
B   721 cactaggcacgactgacgctgaggtacgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaac 800
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A   790 cgatgattactggttgtttgcaatacaccgcaagcgactgagcgaaagcattaagtaatccacctggggagtacgtcggc 869
Q   151 -------------------------------------------------------------------------------- 150
B   801 gatggtaactaggtgtgtgcgacacagagtgcgcgcccaagcgaaagcgataagttacccacctggggagtacgctcgca 880
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A   870 aacgatgaaactcaaaggaattgacgggggcccgcacaagcggtggaacatgtggtttaattcgatgatacgcgaggaac 949
Q   151 -------------------------------------------------------------------------------- 150
B   881 agagtgaaactcaaaggaattgacgggggtccgcacaagcggtggagcatgtggtttaattcgatgatacgcgaggaacc 960
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A   950 cttacctgggtttaaatgggaagtgacaggggtagaaatacctttttcttcggacacttttcaaggtgctgcatggttgt 1029
Q   151 -------------------------------------------------------------------------------- 150
B   961 ttacctgggctcgaatggcctatgacaggcccagagatgggcccttcctcggacataggtcaaggtgctgcatggctgtc 1040
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1030 cgtcagctcgtgccgtgaggtgtcgggttaagtcccataacgagcgcaacccctgttgttagttaccagcatgtaaagat 1109
Q   151 -------------------------------------------------------------------------------- 150
B  1041 gtcagctcgtgccgtgaggtgttgggttaagtcccgcaacgagcgcaacccttgcccctagttgccatcaggtaaagctg 1120
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1110 ggggactctaacaagactgccggtgtaaaccgcgaggaaggtggggatgacgtcaaatcagcacggcccttacatccagg 1189
Q   151 -------------------------------------------------------------------------------- 150
B  1121 gggactctagggggactgcctgcgcaagcagagaggaaggaggggacgatgtcaagtcatcatggcccttacgcccaggg 1200
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1190 gctacacacgtgttacaatggcaggtacaaagggcagctacacagcgatgtgatgctaatctcgaaaacctgtcccagtt 1269
Q   151 -------------------------------------------------------------------------------- 150
B  1201 ctacacacgtgctacaatggcgcatacagagggtagccacctggcgacagggcgccaatctcaaaaagtgcgtctcagtt 1280
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1270 cggattgaagtctgcaacccgacttcatgaagctggaatcgctagtaatcgcgcatcagccatggcgcggtgaatacgtt 1349
Q   151 -------------------------------------------------------------------------------- 150
B  1281 cggatcggggcctgcaactcggccccgtgaagtcggaatcgctagtaatcgcagatcagccatgctgcggtgaatacgtt 1360
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1350 cccgggccttgtacactccgcccgtcaagccatggaagccgggagtacctgaag-------------------------- 1403
Q   151 -------------------------------------------------------------------------------- 150
B  1361 cccgggccttgtacacaccgcccgtcaagccatggaagccgggggcacctgaagtcgggggtaacaacccgcctagggtg 1440
Diffs                                                                                   
Votes                                                                                   
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

A  1404 ---------------------------------------- 1403
Q   151 ---------------------------------------- 150
B  1441 aaactggtaactggggctaagtcgtaacaaggtaaccgta 1480
Diffs                                           
Votes                                           
Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

Ids.  QA 95.3%, QB 91.2%, AB 86.5%, QModel 100.0%, Div. +5.0%
Diffs Left 13: N 0, A 0, Y 13 (100.0%); Right 7: N 0, A 0, Y 7 (100.0%), Score 0.7254
"""

# UCHIME tabular output for single chimeric sequence
single_chimera_ref_tab = """0.7254\t251;size=2;\t4370324\t646991\t4370324\t100.0\t95.3\t91.2\t86.5\t95.3\t13\t0\t0\t7\t0\t0\t4.7\tY
"""


if __name__ == '__main__':
    main()
