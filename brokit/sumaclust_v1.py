#!/usr/bin/env python
"""
Application controller for SumaClust version 1.0
================================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import split, splitext, isdir, exists
from shutil import copy, move
from skbio.app.util import CommandLineApplication, ResultPath
from skbio.app.parameters import ValuedParameter, FlagParameter, FilePath
from skbio.parse.sequences import parse_fasta

class Sumaclust(CommandLineApplication):
    """ SumaClust generic application controller for de novo OTU picking    
    """

    _command = 'sumaclust'
    _command_delimiter = ' '
    _parameters = {
    	# Reference sequence length is the shortest
    	'-l': FlagParameter('-', Name='l', Value=True),
    	# Filepath of the OTU-map
    	'-O': ValuedParameter('-', Name='O', Delimiter=' ', Value=None, IsPath=True),
    	# Number of threads
    	'-p': ValuedParameter('-', Name='p', Delimiter=' ', Value=1, IsPath=False),
    	# Assign sequence to the best matching cluster seed, rather than the first
    	# matching cluster (having >= similarity threshold)
    	'-e': FlagParameter('-', Name='e', Value=False),
    	# Similarity threshold
    	'-t': ValuedParameter('-', Name='t', Delimiter=' ', Value=0.97, IsPath=False)
    }

    _synonyms = {}
    _input_handler = '_input_as_string'
    _supress_stdout = False
    _supress_stderr = False

    def _get_result_paths(self,data):
    	""" Set the result paths
    	"""
    	# directory where all sumaclust's output goes
    	wd = self.WorkingDir
    	result = {}
    	# OTU map (mandatory output)
    	result['OtuMap'] = ResultPath(Path=wd+'sumaclust_otumap.txt',IsWritten=True)
    	# FASTA file of OTUs (mandatory output)
    	result['FastaOTUs'] = ResultPath(Path=wd+'sumaclust_otus.fna',IsWritten=True)
    	return result

	def getHelp(self):
		""" Method that points to documentation
		"""
		help_str =\
		"""
		SumaClust is hosted at:
		http://metabarcoding.org/sumatra/wiki/download

		The following paper should be cited if this resource is used:

		"""

		return help_str


def sumaclust_get_otu_seeds(output_dir=None,seq_path=None):
	"""	Function  : searches for cluster seeds in SumaClust's output FASTA file
	Parameters: output_dir, directory where to store the output file
				seq_path, path to SumaClust's output FASTA file
	Return    : None
	"""

	# Check input files were provided
	if not isdir(output_dir):
		raise ValueError("Output directory to store the FASTA file "
						"of cluster seeds is mandatory input.")
	if not exists(seq_path):
		raise ValueError("The file %s does not exist." % seq_path)

	# Open output file
	extension = splitext(seq_path)[1]
	cluster_centers = output_dir + 'sumaclust_cluster_centers' + extension
	f_out = open(cluster_centers, 'w')

	# Output only the cluster seeds to FASTA file
	f_fasta = open(seq_path, 'U')
	for label, seq in parse_fasta(f_fasta):
		split_label = label.split()
		if "cluster_center=True;" in split_label:
			f_out.write('>'+label+'\n'+seq+'\n')

	f_out.close()
	f_fasta.close()

	return None


def sumaclust_denovo_cluster(seq_path=None,
							output_dir=None,
							shortest_len=True,
							similarity=None,
							threads=1,
							exact=None,
							HALT_EXEC=False
							):
	""" Function  : launch sumaclust de novo OTU picker
		Parameters: 
		Return    : a dictionary of all output files set in 
		            _get_result_paths() with the file 
		            descriptors pointing to an open file
		            as the values
	"""

	# add a '/' to output directory if not already there
	if not output_dir.endswith('/'):
		output_dir = output_dir + '/'

	# Instantiate the object
	sumaclust = Sumaclust(WorkingDir=output_dir, HALT_EXEC=HALT_EXEC)

	# Set the OTU-map filepath
	sumaclust.Parameters['-O'].on(output_dir + 'sumaclust_otumap.txt')

	# Set the similarity threshold
	if similarity is not None:
		sumaclust.Parameters['-t'].on(similarity)

	# Set the option to perform exact clustering (default: False)
	if exact is True:
		sumaclust.Parameters['-e'].on()

	# Launch sumaclust
	# Set the data string to include the read filepath
	# (to be passed as final arguments in the sumaclust command)
	app_result = sumaclust(seq_path)

	# The output OTUs were written to stdout, copy that file
	# to output_dir (temporary solution)
	tmpfile_src = app_result['StdOut'].name 
	tmpfile_dst = output_dir + "sumaclust_otus.fna"

	copy(tmpfile_src,tmpfile_dst)

	# generate a FASTA file of only cluster seeds
	sumaclust_get_otu_seeds(output_dir=output_dir,
	 							seq_path=tmpfile_dst)

	# Return all output files 
	return app_result


























