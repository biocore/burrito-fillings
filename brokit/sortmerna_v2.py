#!/usr/bin/env python
r"""
Application controller for SortMeRNA version 2.0
================================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from os.path import split, splitext
from skbio.app.util import CommandLineApplication, ResultPath
from skbio.app.parameters import ValuedParameter, FlagParameter, FilePath
from glob import glob

class IndexDB(CommandLineApplication):
    """ SortMeRNA generic application controller for building databases    
    """ 
    _command = 'indexdb_rna'
    _command_delimiter = ' '
    _parameters = {
        # Fasta reference file followed by indexed reference
        '--ref': ValuedParameter('--', Name='ref', Delimiter=' ', IsPath=True),

        # Maximum number of positions to store for each unique seed
        '--max_pos': ValuedParameter('--', Name='max_pos', Delimiter=' ', IsPath=False, Value="10000")
    } 

    def _get_result_paths(self, data):
        """ Build the dict of result filepaths
        """
        # access data through self.Parameters so we know it's been cast
        # to a FilePath
        wd = self.WorkingDir
        db_name = (self.Parameters['--ref'].Value).split('/')[-1]

        result = {}
        extensions = ['bursttrie', 'kmer', 'pos', 'stats']
        for extension in extensions:
            for file_path in glob(wd + (db_name + '.' + extension + '*')):
                # this will match e.g. nr.bursttrie_0.dat, nr.bursttrie_1.dat and nr.stats
                key = file_path.split(db_name + '.')[1]
                result_path = ResultPath(Path=file_path, IsWritten=True)
                result[key] = result_path
        return result


def build_database_sortmerna(fasta_path, 
                            max_pos=None, 
                            output_dir=None, 
                            HALT_EXEC=False):
    """ Build sortmerna db from fasta_path; return db name and list of files created

        **If using to create temporary blast databases, you can call
        cogent.util.misc.remove_files(db_filepaths) to clean up all the
        files created by indexdb_rna when you're done with the database.

        fasta_path: path to fasta file of sequences to build database from
        max_pos: maximum positions to store per seed in index 
        output_dir: directory where output should be written
         (default: directory containing fasta_path)
        HALT_EXEC: halt just before running the indexdb_rna command and
         print the command -- useful for debugging
    """

    fasta_dir, fasta_filename = split(fasta_path)
    if not output_dir:
        output_dir = fasta_dir or '.'
        # Will cd to this directory, so just pass the filename
        # so the app is not confused by relative paths
        fasta_path = fasta_filename

    index_basename = splitext(fasta_filename)[0]

    if not output_dir.endswith('/'):
        db_name = output_dir + '/' + index_basename
    else:
        db_name = output_dir + index_basename

    # Instantiate the object
    sdb = IndexDB(WorkingDir=output_dir, HALT_EXEC=HALT_EXEC)
    # The parameter --ref STRING must follow the format where STRING = /path/to/ref.fasta,/path/to/ref.idx
    sdb.Parameters['--ref'].on(fasta_path + ',' + db_name)
    # Override --max_pos parameter default value
    sdb.Parameters['--max_pos'].on(max_pos)

    # Run indexdb_rna
    app_result = sdb()

    db_filepaths = []
    for v in app_result.values():
        try:
            db_filepaths.append(v.name)
        except AttributeError:
            # not a file object, so no path to return
            pass
    return db_name, db_filepaths


class Sortmerna(CommandLineApplication):
    """ SortMeRNA generic application controller for OTU picking    
    """

    _command = 'sortmerna'
    _command_delimiter = ' '
    _parameters = {
        # Fasta or Fastq input reads file
        '--reads': ValuedParameter('--', Name='reads', Delimiter=' ', IsPath=True, Value=None),
    
        # Fasta reference file followed by indexed reference
        '--ref': ValuedParameter('--', Name='ref', Delimiter=' ', IsPath=True, Value=None),
            
        # File path + base name for output files for aligned reads
        '--aligned': ValuedParameter('--', Name='aligned', Delimiter=' ', IsPath=True, Value=None),
            
        # File path + base name for output files for non-aligned reads
        '--other': ValuedParameter('--', Name='other', Delimiter=' ', IsPath=True, Value=None),
            
        # Output Fasta or Fastq file of aligned reads (flag)
        '--fastx': FlagParameter('--', Name='fastx', Value=True),
            
        # Output BLAST alignment file (flag)
        '--blast': ValuedParameter('--', Name='blast', Delimiter=' ', IsPath=False, Value=False),
                        
        # Report the first INT number of alignments
        '--num_alignments': ValuedParameter('--', Name='num_alignments', Delimiter=' ',IsPath=False),
            
        # Report the best INT number of alignments
        '--best': ValuedParameter('--', Name='best', Delimiter=' ',IsPath=False),
            
        # Control the number of reads searched before selecting the best ones
        '--min_lis': ValuedParameter('--', Name='min_lis', Delimiter=' ',IsPath=False),
            
        # Output all reads to SAM and BLAST alignments, a null string for non-aligned reads
        '--print_all_reads': FlagParameter('--',Name='print_all_reads', Value=False),
                                                
        # Number of threads
        '-a': ValuedParameter('-', Name='a', Delimiter=' ', IsPath=False, Value=None),
            
        # E-value threshold
        '-e': ValuedParameter('-', Name='e', Delimiter=' ', IsPath=False, Value=None),
                        
        # Similarity threshold
        '--id': ValuedParameter('--', Name='id', Delimiter=' ',IsPath=False, Value="0.97"),
            
        # Query coverage threshold
        '--coverage': ValuedParameter('--', Name='coverage', Delimiter=' ', IsPath=False, Value="0.97"),
            
        # Output Fasta/Fastq file with reads failing to pass the --id and --coverage thresholds for de novo clustering
        '--de_novo_otu': FlagParameter('--', Name='de_novo_otu', Value=True),
            
        # Output an OTU map
        '--otu_map': FlagParameter('--', Name='otu_map', Value=True)
    }
    _synonyms = {}
    _input_handler = '_input_as_string'
    _supress_stdout = False
    _supress_stderr = False
    _working_dir = None


    def _get_result_paths(self, data):
        ''' Set the result output path
        '''
    
        result = {}
    
        result['OutputAligned'] = ResultPath(Path=self.Parameters['--aligned'.Value],
                                      IsWritten=self.Parameters[--aligned].isOn())
    
        result['OutputNonAligned'] = ResultPath(Path=self.Parameters['--other'.Value],
                                            IsWritten=self.Parameters[--other].isOn())
    
        return result
        
    #def _accept_exit_status(self,exit_status):


# Start reference clustering functions

#def sortmerna_ref_cluster():

