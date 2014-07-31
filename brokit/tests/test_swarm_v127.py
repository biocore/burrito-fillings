#!/usr/bin/env python
"""
Unit tests for the Swarm version 1.2.7 Application controller
=============================================================
"""


from unittest import TestCase, main
import filecmp
from tempfile import mkstemp, mkdtemp
from os import close
from os.path import exists, getsize, join
from shutil import rmtree

from skbio.util.misc import remove_files

from brokit.swarm_v127 import swarm_denovo_cluster


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class SwarmTests(TestCase):
    """ Tests for Swarm version 1.2.7 functionality """

    def setUp(self):
        self.output_dir = mkdtemp()
        self.read_seqs = reads_seqs
        self.reads_seqs_no_derep = reads_seqs_no_derep
        self.derep_reads = derep_reads

        # create temporary file with read sequences defined in read_seqs
        f, self.file_read_seqs = mkstemp(prefix='temp_reads_',
                                         suffix='.fasta')
        close(f)

        # write read sequences to tmp file
        with open(self.file_read_seqs, 'w') as tmp:
            tmp.write(self.read_seqs)

        # create temporary file with non-dereplicated sequences
        f, self.file_reads_seqs_no_derep = mkstemp(
            prefix='temp_reads_no_derep',
            suffix='.txt')
        close(f)

        # write non-dereplicated sequences to file
        with open(self.file_reads_seqs_no_derep, 'w') as tmp:
            tmp.write(self.reads_seqs_no_derep)

        # create temporary file with de-replicated sequences
        f, self.file_derep_reads = mkstemp(prefix='temp_derep_reads_',
                                           suffix='.fasta')
        close(f)

        # write de-replicated sequences to file
        with open(self.file_derep_reads, 'w') as tmp:
            tmp.write(self.derep_reads)

        # list of files to remove
        self.files_to_remove = [self.file_read_seqs,
                                self.file_reads_seqs_no_derep,
                                self.file_derep_reads]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.output_dir)

    def test_default_param(self):
        """ Swarm should return the correct clusters using
            default inputs
        """
        clusters = swarm_denovo_cluster(seq_path=self.file_read_seqs,
                                        output_dir=self.output_dir,
                                        d=1,
                                        threads=1,
                                        prefilter_identical_sequences=True)

        # Check the returned clusters list of lists is as expected
        expected_clusters = [['s1_630', 's1_4572', 's1_5748',
                              's1_13961', 's1_8744', 's1_8623',
                              's1_7634', 's1_6846', 's1_3750',
                              's1_2369'],
                             ['s1_8977', 's1_10439', 's1_12366',
                              's1_15985', 's1_21935', 's1_8592',
                              's1_4677', 's1_14735', 's1_11650',
                              's1_11001'],
                             ['s1_844', 's1_1886', 's1_5347',
                              's1_5737', 's1_7014', 's1_7881',
                              's1_8615', 's1_7040', 's1_6200',
                              's1_1271']]

        # Should be 3 clusters
        self.assertEqual(len(clusters), 3)

        # List of actual clusters matches list of expected clusters
        for actual_cluster, expected_cluster in zip(clusters,
                                                    expected_clusters):
            actual_cluster.sort()
            expected_cluster.sort()
            self.assertEqual(actual_cluster, expected_cluster)

    def test_prefilter_identical_sequences_false(self):
        """ The option prefilter_identical_sequences can be
            set to false only if the input FASTA file has been
            de-replicated to include the abundance information
            after the last underscore '_' in the sequence label,
            otherwise a ValueError is thrown.
            Note: "prefilter_identical_sequences=False" should
            be used with caution since there's no easy way for
            us to verify that the number after the final '_' in
            the label is indeed the abundance for dereplicated
            input sequences (rather than part of the sequence ID),
            we leave this to the user to check
        """

        self.assertRaises(ValueError,
                          swarm_denovo_cluster,
                          seq_path=self.file_reads_seqs_no_derep,
                          output_dir=self.output_dir,
                          d=1,
                          threads=1,
                          prefilter_identical_sequences=False)

    def test_pass_dereplicated_reads(self):
        """ Reads will not be de-replicated if
            'prefilter_identical_sequences=False' and
            the sequence labels appear to have abundance
            as the number after the final '_'.
            Note: "prefilter_identical_sequences=False" should
            be used with caution since there's no easy way for
            us to verify that the number after the final '_' in
            the label is indeed the abundance for dereplicated
            input sequences (rather than part of the sequence ID),
            we leave this to the user to check
        """
        clusters = swarm_denovo_cluster(seq_path=self.file_derep_reads,
                                        output_dir=self.output_dir,
                                        d=1,
                                        threads=1,
                                        prefilter_identical_sequences=False)

        # Check the returned clusters list of lists is as expected
        expected_clusters = [['QiimeExactMatch.s1_630',
                              'QiimeExactMatch.s1_8744',
                              'QiimeExactMatch.s1_8623',
                              'QiimeExactMatch.s1_7634',
                              'QiimeExactMatch.s1_6846',
                              'QiimeExactMatch.s1_3750',
                              'QiimeExactMatch.s1_2369'],
                             ['QiimeExactMatch.s1_8977',
                              'QiimeExactMatch.s1_8592',
                              'QiimeExactMatch.s1_4677',
                              'QiimeExactMatch.s1_14735',
                              'QiimeExactMatch.s1_11650',
                              'QiimeExactMatch.s1_11001'],
                             ['QiimeExactMatch.s1_844',
                              'QiimeExactMatch.s1_8615',
                              'QiimeExactMatch.s1_7040',
                              'QiimeExactMatch.s1_6200',
                              'QiimeExactMatch.s1_1271']]

        # Should be 3 clusters
        self.assertEqual(len(clusters), 3)

        # List of actual clusters matches list of expected clusters
        for actual_cluster, expected_cluster in zip(clusters,
                                                    expected_clusters):
            actual_cluster.sort()
            expected_cluster.sort()
            self.assertEqual(actual_cluster, expected_cluster)

    def test_seq_path(self):
        """ Swarm should raise a ValueError if the sequences
            filepath was not provided
        """

        self.assertRaises(ValueError,
                          swarm_denovo_cluster,
                          seq_path=None,
                          output_dir=self.output_dir,
                          d=1,
                          threads=1,
                          prefilter_identical_sequences=True)

    def test_output_dir(self):
        """ Swarm should raise a ValueError if the output
            directory was not provided
        """

        self.assertRaises(ValueError,
                          swarm_denovo_cluster,
                          seq_path=self.file_read_seqs,
                          output_dir=None,
                          d=1,
                          threads=1,
                          prefilter_identical_sequences=True)

    def test_negative_resolution(self):
        """ Swarm should raise a ValueError if the resolution
            is negative
        """

        self.assertRaises(ValueError,
                          swarm_denovo_cluster,
                          seq_path=self.file_read_seqs,
                          output_dir=self.output_dir,
                          d=-2,
                          threads=1,
                          prefilter_identical_sequences=True)

    def test_negative_threads(self):
        """ Swarm should raise a ValueError if number of threads
            is negative
        """

        self.assertRaises(ValueError,
                          swarm_denovo_cluster,
                          seq_path=self.file_read_seqs,
                          output_dir=self.output_dir,
                          d=1,
                          threads=-2,
                          prefilter_identical_sequences=True)

# Reads to cluster
# there are 30 reads representing 3 species (gives 3 clusters)
reads_seqs = """>s1_630 reference=1049393 amplicon=complement(497..788)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_2369 reference=1049393 amplicon=complement(497..788) errors=73%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTAGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_3750 reference=1049393 amplicon=complement(497..788) errors=100%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCA
>s1_4572 reference=1049393 amplicon=complement(497..788)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_5748 reference=1049393 amplicon=complement(497..788)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_6846 reference=1049393 amplicon=complement(497..788) errors=67%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCATAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_7634 reference=1049393 amplicon=complement(497..788) errors=99%T
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTTG
>s1_8623 reference=1049393 amplicon=complement(497..788) errors=17-
GTGCCAGCAGCCGCGGAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_8744 reference=1049393 amplicon=complement(497..788) errors=62%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGAGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_13961 reference=1049393 amplicon=complement(497..788)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_4677 reference=4382408 amplicon=complement(487..778) errors=74%T
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGTGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_8592 reference=4382408 amplicon=complement(487..778) errors=95+A
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAAGCCCA
>s1_8977 reference=4382408 amplicon=complement(487..778)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_10439 reference=4382408 amplicon=complement(487..778)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_11001 reference=4382408 amplicon=complement(487..778) errors=91%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGGGAAAGCCCA
>s1_11650 reference=4382408 amplicon=complement(487..778) errors=78-
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCGTAAGTCAGAGGTGAAAGCCCA
>s1_12366 reference=4382408 amplicon=complement(487..778)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_14735 reference=4382408 amplicon=complement(487..778) errors=94%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGACAGCCCA
>s1_15985 reference=4382408 amplicon=complement(487..778)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_21935 reference=4382408 amplicon=complement(487..778)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_844 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_1271 reference=129416 amplicon=complement(522..813) errors=94%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGACAGCCCA
>s1_1886 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_5347 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_5737 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_6200 reference=129416 amplicon=complement(522..813) errors=92%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTCAAAGCCCA
>s1_7014 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_7040 reference=129416 amplicon=complement(522..813) errors=40%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAGTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_7881 reference=129416 amplicon=complement(522..813)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_8615 reference=129416 amplicon=complement(522..813) errors=81%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTGAGTCAGATGTGAAAGCCCA
"""

# not de-replicated sequences
reads_seqs_no_derep = """>s1 reference=1049393 amplicon=complement(497..788)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s2 reference=1049393 amplicon=complement(497..788) errors=73%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTAGCGGGGTAAGTCAGGTGTGAAATCTCG
>s3 reference=1049393 amplicon=complement(497..788) errors=100%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCA
>s4 reference=1049393 amplicon=complement(497..788)
"""

# de-replicated sequences
derep_reads = """>QiimeExactMatch.s1_630_4
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>QiimeExactMatch.s1_2369_1
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTAGCGGGGTAAGTCAGGTGTGAAATCTCG
>QiimeExactMatch.s1_3750_1
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCA
>QiimeExactMatch.s1_6846_1
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCATAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>QiimeExactMatch.s1_7634_1
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTTG
>QiimeExactMatch.s1_8623_1
GTGCCAGCAGCCGCGGAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>QiimeExactMatch.s1_8744_1
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGAGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>QiimeExactMatch.s1_4677_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGTGTCTGTAAGTCAGAGGTGAAAGCCCA
>QiimeExactMatch.s1_8592_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAAGCCCA
>QiimeExactMatch.s1_8977_5
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>QiimeExactMatch.s1_11001_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGGGAAAGCCCA
>QiimeExactMatch.s1_11650_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCGTAAGTCAGAGGTGAAAGCCCA
>QiimeExactMatch.s1_14735_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGACAGCCCA
>QiimeExactMatch.s1_844_6
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>QiimeExactMatch.s1_1271_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGACAGCCCA
>QiimeExactMatch.s1_6200_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTCAAAGCCCA
>QiimeExactMatch.s1_7040_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAGTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>QiimeExactMatch.s1_8615_1
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTGAGTCAGATGTGAAAGCCCA
"""

if __name__ == '__main__':
    main()
