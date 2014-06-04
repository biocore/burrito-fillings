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


from string import strip
from skbio.app.util import (CommandLineApplication, ResultPath,
                            get_tmp_filename, guess_input_handler,
                            ApplicationNotFoundError)


class Sortmerna(CommandLineApplication):
    """ SortMeRNA generic application controller
    
    """

    _command = 'sortmerna'
    _command_delimiter = ' '
    _parameters = {

        # Fasta or Fastq input reads file
        '--reads': ValuedParameter('--',Name='reads', Delimiter=' ',IsPath=True),
    
        # Fasta reference file followed by indexed reference
        '--ref': ValuedParameter('--',Name='reference', Delimiter=' ',IsPath=True),
            
        # File path + base name for output files for aligned reads
        '--aligned': ValuedParameter('--',Name='aligned', Delimiter=' ',IsPath=True),
            
        # File path + base name for output files for non-aligned reads
        '--other': ValuedParameter('--',Name='other', Delimiter=' ',IsPath=True),
            
        # Output Fasta or Fastq file of aligned reads (flag)
        '--fastx': FlagParameter('--', Name='fastx'),
            
        # Output SAM alignment file (flag)
        '--sam': FlagParameter('--', Name='sam'),
            
        # For SAM alignment, output SQ tags (flag)
        '--SQ': FlagParameter('--', Name='SQ'),
            
        # Output BLAST alignment file (flag)
        '--blast': ValuedParameter('--', Name='blast',Delimiter=' ',IsPath=False),
            
        # Output log file with all alignment statistics (flag)
        '--log': FlagParameter('--', Name='log'),
            
        # Report the first INT number of alignments
        '--num_alignments': ValuedParameter('--',Name='num_alignments', Delimiter=' ',IsPath=False),
            
        # Report the best INT number of alignments
        '--best': ValuedParameter('--',Name='best', Delimiter=' ',IsPath=False),
            
        # Control the number of reads searched before selecting the best ones
        '--min_lis': ValuedParameter('--',Name='min_lis', Delimiter=' ',IsPath=False),
            
        # Output all reads to SAM and BLAST alignments, a null string for non-aligned reads
        '--print_all_reads': FlagParameter('--',Name='print_all_reads'),
            
        # Output paired reads together in aligned Fasta/Fastq file (even if one of them didn't align)
        '--paired_in': FlagParameter('--',Name='paired_in'),
                
        # Output paired reads together in non-aligned Fasta/Fastq file (even if one of them aligned)
        '--paired_out': FlagParameter('--',Name='paired_out'),
            
        # Smith-Waterman match score (positive integer)
        '--match': ValuedParameter('--',Name='match', Delimiter=' ',IsPath=False),
            
        # Smith-Waterman mismatch penalty (negative integer)
        '--mismatch': ValuedParameter('--',Name='mismatch', Delimiter=' ',IsPath=False),
            
        # Smith-Waterman gap open penalty (positive integer)
        '--gap_open': ValuedParameter('--',Name='gap_open', Delimiter=' ',IsPath=False),
            
        # Smith-Waterman gap extend penalty (positive integer)
        '--gap_ext': ValuedParameter('--',Name='gap_ext', Delimiter=' ',IsPath=False),
            
        # Smith-Waterman penalty for ambigious N character
        '-N': ValuedParameter('-',Name='N', Delimiter=' ',IsPath=False),
            
        # Search only the forward strand
        '-F': ValuedParameter('-',Name='F', Delimiter=' ',IsPath=False),
            
        # Search only the reverse strand
        '-R': ValuedParameter('-',Name='R', Delimiter=' ',IsPath=False),
            
        # Number of threads
        '-a': ValuedParameter('-',Name='a', Delimiter=' ',IsPath=False),
            
        # E-value threshold
        '-e': ValuedParameter('-',Name='e', Delimiter=' ',IsPath=False),
            
        # Number of Mbytes for loading the reads into memory
        '-m': ValuedParameter('-',Name='m', Delimiter=' ',IsPath=False),
            
        # Verbose
        '-v': FlagParameter('-',Name='v'),
            
        # Similarity threshold
        '--id': ValuedParameter('--',Name='id', Delimiter=' ',IsPath=False),
            
        # Query coverage threshold
        '--coverage': ValuedParameter('--',Name='coverage', Delimiter=' ',IsPath=False),
            
        # Output Fasta/Fastq file with reads failing to pass the --id and --coverage thresholds for de novo clustering
        '--de_novo_otu': FlagParameter('--',Name='de_novo_otu'),
            
        # Output an OTU map
        '--otu_map': FlagParameter('--',Name='otu_map'),
            
        # Intervals at which to place the seed on the read
        '--passes': ValuedParameter('--',Name='passes', Delimiter=' ',IsPath=False),
            
        # number (or %) of nucleotides on each end of the read to consider during alignment
        '--edges': ValuedParameter('--',Name='edges', Delimiter=' ',IsPath=False),
            
        # Number of seeds to match before searching for candidate LIS
        '--num_seeds': ValuedParameter('--',Name='num_seeds', Delimiter=' ',IsPath=False),
            
        # Search for all 0-error and 1-error seeds rather than stopping after finding 0-error match
        '--full_search': FlagParameter('--',Name='full_search'),
            
        # Attach process ID to the end of output files
        '--pid': FlagParameter('--',Name='pid'),
            
        # Print help page
        '-h': FlagParameter('-',Name='pid'),
            
        # Print Sortmerna version
        '--version': FlagParameter('--',Name='version'),

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
    
    def _error_on_missing_application(self,params):
        """ Raise an ApplicationNotFoundError if sortmerna or indexdb_rna are not accessible
        """
    
        if not app_path('sortmerna'):
            raise ApplicationNotFoundError, "Cannot find sortmerna. Is it installed and in your path?"
        
        elif not app_path('indexdb_rna'):
            raise ApplicationNotFoundError, "Cannot find indexdb_rna. Is it installed and in your path?"
    
    def _accept_exit_status(self,exit_status):


# Start reference clustering functions

def sortmerna_v2_build_index(
):

def sortmerna_v2_ref_cluster(
):
