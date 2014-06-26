#!/usr/bin/env python
"""
Unit tests for the SumaClust version 1.0 Application controller
===============================================================
"""


from unittest import TestCase, main
import filecmp
from tempfile import mkstemp, mkdtemp
from os import close
from os.path import exists, getsize, join
from shutil import rmtree

from skbio.util.misc import remove_files

from brokit.sumaclust_v1 import (sumaclust_denovo_cluster)


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class SumaclustV1Tests(TestCase):
    """ Tests for Sumaclust version 2.0 functionality """

    def setUp(self):

        self.output_dir = mkdtemp()
        self.read_seqs = reads_seqs
        self.expected_otumap = expected_otumap

        # create temporary file with read sequences defined in read_seqs
        f, self.file_read_seqs = mkstemp(prefix='temp_reads_',
                                         suffix='.fasta')
        close(f)

        # write read sequences to tmp file
        with open(self.file_read_seqs, 'w') as tmp:
            tmp.write(self.read_seqs)

        # create temporary file with final OTU map
        f, self.file_otumap = mkstemp(prefix='temp_otumap',
                                      suffix='.txt')
        close(f)

        # write OTU map to tmp file
        with open(self.file_otumap, 'w') as tmp:
            tmp.write(self.expected_otumap)

        # list of files to remove
        self.files_to_remove = [self.file_read_seqs, self.file_otumap]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.output_dir)

    def check_clusters(self,
                       clusters=None,
                       result_path=None):

        # Check the OTU map was output with the correct size
        self.assertTrue(exists(result_path))

        # Checkout output file has the correct size
        size = getsize(result_path)
        self.assertTrue(size, 270)

        # Check file identical to expected result (the hash should
        # be identical)
        self.assertTrue(filecmp.cmp(result_path, self.file_otumap))

        # Check the returned clusters list of lists is as expected
        expected_clusters = [['s1_844', 's1_1886', 's1_5347', 's1_5737',
                              's1_7014', 's1_7881', 's1_7040', 's1_6200',
                              's1_1271', 's1_8615'],
                             ['s1_8977', 's1_10439', 's1_12366', 's1_15985',
                              's1_21935', 's1_11650', 's1_11001', 's1_8592',
                              's1_14735', 's1_4677'],
                             ['s1_630', 's1_4572', 's1_5748', 's1_13961',
                              's1_2369', 's1_3750', 's1_7634', 's1_8623',
                              's1_8744', 's1_6846']]

        # Should be 3 clusters
        self.assertEqual(len(clusters), 3)

        # List of actual clusters matches list of expected clusters
        for actual_cluster, expected_cluster in zip(clusters,
                                                    expected_clusters):
            actual_cluster.sort()
            expected_cluster.sort()
            self.assertEqual(actual_cluster, expected_cluster)

    def test_empty_seq_path(self):
        """ SumaClust should return a ValueError
            if empty sequence path is passed
        """
        result_path = join(self.output_dir, "sumaclust_otus.txt")

        self.assertRaises(ValueError,
                          sumaclust_denovo_cluster,
                          seq_path=None,
                          result_path=result_path)

    def test_empty_result_path(self):
        """ SumaClust should return a ValueError
            if empty result path is passed
        """
        self.assertRaises(ValueError,
                          sumaclust_denovo_cluster,
                          seq_path=self.file_read_seqs,
                          result_path=None)

    def test_negative_threads(self):
        """ SumaClust should raise ValueError
            on negative number of threads
        """
        result_path = join(self.output_dir, "sumaclust_otus.txt")

        self.assertRaises(ValueError,
                          sumaclust_denovo_cluster,
                          seq_path=self.file_read_seqs,
                          result_path=result_path,
                          shortest_len=True,
                          similarity=0.97,
                          threads=-2)

    def test_positive_threads(self):
        """ SumaClust's actual clusters should match
            the exact clusters when using multithreading
        """
        result_path = join(self.output_dir, "sumaclust_otus_exact.txt")
        clusters = sumaclust_denovo_cluster(seq_path=self.file_read_seqs,
                                            result_path=result_path,
                                            shortest_len=True,
                                            similarity=0.97,
                                            threads=3,
                                            exact=True)

        self.files_to_remove.append(result_path)

        self.check_clusters(clusters, result_path)

    def test_exact_clustering(self):
        """ SumaClust's actual clusters should match
            the exact clusters when using the exact option
        """
        result_path = join(self.output_dir, "sumaclust_otus_exact.txt")
        clusters = sumaclust_denovo_cluster(seq_path=self.file_read_seqs,
                                            result_path=result_path,
                                            shortest_len=True,
                                            similarity=0.97,
                                            threads=1,
                                            exact=True)

        self.files_to_remove.append(result_path)

        self.check_clusters(clusters, result_path)

    def test_shortest_len_clustering(self):
        """ SumaClust's actual clusters should match
            the exact clusters when not using the
            shortest len option
        """
        result_path = join(self.output_dir, "sumaclust_otus_exact.txt")
        clusters = sumaclust_denovo_cluster(seq_path=self.file_read_seqs,
                                            result_path=result_path,
                                            shortest_len=False,
                                            similarity=0.97,
                                            threads=1,
                                            exact=True)

        self.files_to_remove.append(result_path)

        self.check_clusters(clusters, result_path)

    def test_sumaclust_denovo_cluster(self):
        """ Test de novo clustering with SumaClust """

        result_path = join(self.output_dir, "sumaclust_otus.txt")

        clusters = sumaclust_denovo_cluster(seq_path=self.file_read_seqs,
                                            result_path=result_path)

        self.files_to_remove.append(result_path)

        self.check_clusters(clusters, result_path)


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

# The expected OTU map from clustering the read_seqs from 3 species
# (using exact and non-exact options)
expected_otumap = """s1_844\ts1_844\ts1_1886\ts1_5347\ts1_5737\ts1_7014\ts1_7881\ts1_7040\ts1_6200\ts1_1271\ts1_8615
s1_8977\ts1_8977\ts1_10439\ts1_12366\ts1_15985\ts1_21935\ts1_11650\ts1_11001\ts1_8592\ts1_14735\ts1_4677
s1_13961\ts1_13961\ts1_4572\ts1_5748\ts1_630\ts1_2369\ts1_3750\ts1_7634\ts1_8623\ts1_8744\ts1_6846
"""

if __name__ == '__main__':
    main()
