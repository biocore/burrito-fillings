#!/usr/bin/env python
"""
Unit tests for the SumaClust version 1.0 Application controller
===============================================================
"""


from unittest import TestCase, main
from brokit.sumaclust_v1 import (Sumaclust,
								 sumaclust_denovo_cluster)
from tempfile import mkstemp, mkdtemp
from skbio.util.misc import remove_files
from os import close
from os.path import abspath, exists, getsize
from skbio.parse.sequences import parse_fasta
import re
import filecmp


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


# Reads to cluster
# there are 30 reads representing 3 species (gives 3 clusters)
reads_seqs = """>s1_630 reference=1049393 amplicon=complement(497..788) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_2369 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=73%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTAGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_3750 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=100%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCA
>s1_4572 reference=1049393 amplicon=complement(497..788) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_5748 reference=1049393 amplicon=complement(497..788) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_6846 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=67%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCATAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_7634 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=99%T
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTTG
>s1_8623 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=17-
GTGCCAGCAGCCGCGGAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_8744 reference=1049393 amplicon=complement(497..788) position=complement(193..292) errors=62%A
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGAGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_13961 reference=1049393 amplicon=complement(497..788) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGCGGGGTAAGTCAGGTGTGAAATCTCG
>s1_4677 reference=4382408 amplicon=complement(487..778) position=complement(193..292) errors=74%T
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGTGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_8592 reference=4382408 amplicon=complement(487..778) position=complement(193..292) errors=95+A
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAAGCCCA
>s1_8977 reference=4382408 amplicon=complement(487..778) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_10439 reference=4382408 amplicon=complement(487..778) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_11001 reference=4382408 amplicon=complement(487..778) position=complement(193..292) errors=91%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGGGAAAGCCCA
>s1_11650 reference=4382408 amplicon=complement(487..778) position=complement(193..292) errors=78-
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCGTAAGTCAGAGGTGAAAGCCCA
>s1_12366 reference=4382408 amplicon=complement(487..778) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_14735 reference=4382408 amplicon=complement(487..778) position=complement(193..292) errors=94%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGACAGCCCA
>s1_15985 reference=4382408 amplicon=complement(487..778) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_21935 reference=4382408 amplicon=complement(487..778) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTCCAAGCGTTGTCCGGAATCACTGGGTGTAAAGGGTGCGTAGGCGGGTCTGTAAGTCAGAGGTGAAAGCCCA
>s1_844 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_1271 reference=129416 amplicon=complement(522..813) position=complement(193..292) errors=94%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGACAGCCCA
>s1_1886 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_5347 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_5737 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_6200 reference=129416 amplicon=complement(522..813) position=complement(193..292) errors=92%C
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTCAAAGCCCA
>s1_7014 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_7040 reference=129416 amplicon=complement(522..813) position=complement(193..292) errors=40%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAGTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_7881 reference=129416 amplicon=complement(522..813) position=complement(193..292)
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTAAGTCAGATGTGAAAGCCCA
>s1_8615 reference=129416 amplicon=complement(522..813) position=complement(193..292) errors=81%G
GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGCTTTGTGAGTCAGATGTGAAAGCCCA
""" 

# The expected OTU map from clustering the read_seqs from 3 species
expected_otumap = """s1_844\ts1_844\ts1_1886\ts1_5347\ts1_5737\ts1_7014\ts1_7881\ts1_7040\ts1_6200\ts1_1271\ts1_8615
s1_8977\ts1_8977\ts1_10439\ts1_12366\ts1_15985\ts1_21935\ts1_11650\ts1_11001\ts1_8592\ts1_14735\ts1_4677
s1_13961\ts1_13961\ts1_4572\ts1_5748\ts1_630\ts1_2369\ts1_3750\ts1_7634\ts1_8623\ts1_8744\ts1_6846
"""

class SumaclustV1Tests(TestCase):
	""" Tests for Sumaclust version 2.0 functionality """

	def setUp(self):

		self.output_dir = '/tmp/'
		self.read_seqs = reads_seqs
		self.expected_otumap = expected_otumap

		# create temporary file with read sequences defined in read_seqs
		f, self.file_read_seqs = mkstemp(prefix='temp_reads_',
										 suffix='.fasta')
		close(f)

		# write read sequences to tmp file
		tmp = open(self.file_read_seqs, 'w')
		tmp.write(self.read_seqs)
		tmp.close()

		# create temporary file with final OTU map
		f, self.file_otumap = mkstemp(prefix='temp_otumap',
										suffix='.txt')
		close(f)

		# write OTU map to tmp file
		tmp = open(self.file_otumap, 'w')
		tmp.write(self.expected_otumap)
		tmp.close()

		# list of files to remove
		self.files_to_remove = [self.file_read_seqs, self.file_otumap]

	def tearDown(self):
		pass
		#remove_files(self.files_to_remove)

	def test_sumaclust_denovo_cluster(self):
		""" Test de novo clustering with SumaClust """

		result_path = self.output_dir + "sumaclust_otus.txt"

		output_files = sumaclust_denovo_cluster(seq_path=self.file_read_seqs,
												result_path=result_path,
												shortest_len=True,
												similarity=None,
												threads=1,
												exact=False,
												HALT_EXEC=False)

		# Check the OTU map was output with the correct size
		otumap = output_files['OtuMap'].name
		self.assertTrue(exists(otumap))

		# Checkout output file has the correct size
		size = getsize(otumap)
		self.assertTrue(size,270)

		# Check file identical to expected result (the hash should be identical)
		self.assertTrue(filecmp.cmp(otumap,self.file_otumap))

if __name__ == '__main__':
    main()

