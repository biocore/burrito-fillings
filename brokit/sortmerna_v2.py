#!/usr/bin/env python
"""
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


from os.path import split, splitext, dirname
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
        '--max_pos': ValuedParameter('--', Name='max_pos', Delimiter=' ', IsPath=False, Value="10000"),

        # tmp folder for storing unique L-mers (prior to calling CMPH in indexdb_rna),
        # this tmp file is removed by indexdb_rna after it is not used any longer
        '--tmpdir': ValuedParameter('--', Name='tmpdir', Delimiter=' ', IsPath=True)
    } 

    def _get_result_paths(self, data):
        """ Build the dict of result filepaths
        """
        wd = self.WorkingDir

        db_name = (self.Parameters['--ref'].Value).split('/')[-1]

        result = {}
        extensions = ['bursttrie', 'kmer', 'pos', 'stats']
        for extension in extensions:
            for file_path in glob(wd + (db_name + '.' + extension + '*')):
                # this will match e.g. nr.bursttrie_0.dat, nr.bursttrie_1.dat and nr.stats
                key = file_path.split(db_name + '.')[1]
                result[key] = ResultPath(Path=file_path, IsWritten=True)                 
        return result



    def _accept_exit_status(self, exit_status):
        """ Test for acceptable exit status
        """
        return exit_status == 0


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

    # Set temporary directory
    sdb.Parameters['--tmpdir'].on(output_dir)

    # Override --max_pos parameter 
    if max_pos is not None:
        sdb.Parameters['--max_pos'].on(max_pos)

    # Run indexdb_rna
    app_result = sdb()

    # Return all output files (by indexdb_rna) as a list,
    # first however remove the StdErr and StdOut filepaths
    # as they files will be destroyed at the exit from
    # this function (IndexDB is a local instance)
    db_filepaths = []
    for key in app_result:
        try:
            if (key is not 'StdErr' and \
                key is not 'StdOut'):
                db_filepaths.append(app_result[key].name)
        except AttributeError:
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

        # Output log file with parameters used to launch sortmerna and statistics on final results
        '--log': FlagParameter('--', Name='log', Value=True),
                        
        # Output Fasta or Fastq file of aligned reads (flag)
        '--fastx': FlagParameter('--', Name='fastx', Value=True),
            
        # Output BLAST alignment file (flag)
        '--blast': ValuedParameter('--', Name='blast', Delimiter=' ', IsPath=False, Value=None),
                                    
        # Report the best INT number of alignments
        '--best': ValuedParameter('--', Name='best', Delimiter=' ',IsPath=False, Value="1"),
                                                            
        # Number of threads
        '-a': ValuedParameter('-', Name='a', Delimiter=' ', IsPath=False, Value="1"),
            
        # E-value threshold
        '-e': ValuedParameter('-', Name='e', Delimiter=' ', IsPath=False, Value="1"),
                        
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


    def _get_result_paths(self,data):
        """ Set the result paths """

        result = {}

        # get the file extension of the reads file (sortmerna internally outputs 
        # all results with this extension)
        fileExtension = splitext(self.Parameters['--reads'].Value)[1]

        # at this point the parameter --aligned should be set as 
        # sortmerna will not run without it
        if not self.Parameters['--aligned'].isOn():
            print "Error: "
            sys.exit(0)

        # file base name for aligned reads
        output_base = self.Parameters['--aligned'].Value

        # Blast alignments
        if self.Parameters['--blast'].isOn():
            result['BlastAlignments'] = ResultPath(Path=output_base + '.blast',
                                                   IsWritten=True)

        # OTU map (mandatory output)
        result['OtuMap'] = ResultPath(Path=output_base + '_otus.txt',
                                      IsWritten=True)

        # FASTA file of sequences in the OTU map (madatory output)
        result['FastaMatches'] = ResultPath(Path=output_base + fileExtension,
                                            IsWritten=True)

        # FASTA file of sequences not in the OTU map (mandatory output)
        result['FastaForDenovo'] = ResultPath(Path=output_base + '_denovo' + fileExtension,
                                              IsWritten=True)

        # Log file
        result['LogFile'] = ResultPath(Path=output_base + '.log',
                                        IsWritten=True)
    
        return result

    def _accept_exit_status(self, exit_status):
        """ Test for acceptable exit status
        """
        return exit_status == 0

    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
            """
        SortMeRNA is hosted at:
        http://bioinfo.lifl.fr/RNA/sortmerna/
        https://github.com/biocore/sortmerna

        The following paper should be cited if this resource is used:

        Kopylova, E., Noe L. and Touzet, H.,
        SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data, 
        Bioinformatics (2012) 28(24)
        """
        return help_str
        



# Start reference clustering functions

def sortmerna_ref_cluster(seq_path=None,
                          sortmerna_db=None,
                          refseqs_fp=None,
                          result_path=None,
                          max_e_value=None,
                          similarity=None,
                          coverage=None,
                          threads=None,
                          tabular=False,
                          best=None,
                          HALT_EXEC=False
                          ):
    ''' Function  : Launch sortmerna OTU picker 
        Parameters: sortmerna_db, refseqs_fp, seq_path and output_dir are mandatory and must already exist 
        Return    : a dictionary of all output files set in _get_result_paths() with the file descriptors
                    pointing to an open file as the values
    '''

    # Instantiate the object
    smr = Sortmerna(HALT_EXEC=HALT_EXEC)

    # Set input reads path
    if seq_path is not None:
        smr.Parameters['--reads'].on(seq_path)
    else:
        print "Error: an read file is mandatory input. " 
        sys.exit(1)

    # Set the input reference sequence + indexed database path
    if sortmerna_db is not None:
        smr.Parameters['--ref'].on(refseqs_fp + ',' + sortmerna_db)
    else:
        print "Error: a indexed database for reference set %s must already exist. " % refseqs_fp
        print "       Use indexdb_rna to index the database."
        sys.exit(1)

    # Set output results path (for Blast alignments, clusters and failures)
    output_dir = dirname(result_path)
    if output_dir is not None:
        if not output_dir.endswith('/'):
            output_file = output_dir + '/sortmerna_otus'
        else:
            output_file = output_dir + 'sortmerna_otus'
        smr.Parameters['--aligned'].on(output_file)

    # Set E-value threshold
    if max_e_value is not None:
        smr.Parameters['-e'].on(max_e_value)

    # Set similarity threshold
    if similarity is not None:
        smr.Parameters['--id'].on(similarity)

    # Set query coverage threshold
    if coverage is not None:
        smr.Parameters['--coverage'].on(coverage)

    # Set number of best alignments to output
    if best is not None:
        smr.Parameters['--best'].on(best)

    if tabular:
        smr.Parameters['--blast'].on("3")

    # Run sortmerna
    app_result = smr()

    return app_result

