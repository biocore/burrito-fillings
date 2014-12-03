#!/usr/bin/env python

"""
Application controller for SortMeRNA version 2.0
================================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from os.path import split, splitext, dirname, join
from os import getpid
from glob import glob
import re

from burrito.util import CommandLineApplication, ResultPath
from burrito.parameters import ValuedParameter, FlagParameter
from skbio.parse.sequences import parse_fasta


class IndexDB(CommandLineApplication):
    """ SortMeRNA generic application controller for building databases
    """
    _command = 'indexdb_rna'
    _command_delimiter = ' '
    _parameters = {
        # Fasta reference file followed by indexed reference
        # (ex. /path/to/refseqs.fasta,/path/to/refseqs.idx)
        '--ref': ValuedParameter('--', Name='ref', Delimiter=' ', IsPath=True),

        # Maximum number of positions to store for each unique seed
        '--max_pos': ValuedParameter('--', Name='max_pos', Delimiter=' ',
                                     IsPath=False, Value="10000"),

        # tmp folder for storing unique L-mers (prior to calling CMPH
        # in indexdb_rna), this tmp file is removed by indexdb_rna
        # after it is not used any longer
        '--tmpdir': ValuedParameter('--', Name='tmpdir', Delimiter=' ',
                                    IsPath=True)
    }

    def _get_result_paths(self, data):
        """ Build the dict of result filepaths
        """
        # get the filepath of the indexed database (after comma)
        # /path/to/refseqs.fasta,/path/to/refseqs.idx
        #                        ^------------------^
        db_name = (self.Parameters['--ref'].Value).split(',')[1]

        result = {}
        extensions = ['bursttrie', 'kmer', 'pos', 'stats']
        for extension in extensions:
            for file_path in glob("%s.%s*" % (db_name, extension)):
                # this will match e.g. nr.bursttrie_0.dat, nr.bursttrie_1.dat
                # and nr.stats
                key = file_path.split(db_name + '.')[1]
                result[key] = ResultPath(Path=file_path, IsWritten=True)
        return result


def build_database_sortmerna(fasta_path,
                             max_pos=None,
                             output_dir=None,
                             HALT_EXEC=False):
    """ Build sortmerna db from fasta_path; return db name
        and list of files created

        Parameters
        ----------
        fasta_path : string
            path to fasta file of sequences to build database.
        max_pos : integer, optional
            maximum positions to store per seed in index
            [default: 10000].
        output_dir : string, optional
            directory where output should be written
            [default: same directory as fasta_path]
        HALT_EXEC : boolean, optional
            halt just before running the indexdb_rna command
            and print the command -- useful for debugging
            [default: False].

        Return
        ------
        db_name : string
            filepath to indexed database.
        db_filepaths : list
            output files by indexdb_rna
    """

    if fasta_path is None:
        raise ValueError("Error: path to fasta reference "
                         "sequences must exist.")

    fasta_dir, fasta_filename = split(fasta_path)
    if not output_dir:
        output_dir = fasta_dir or '.'
        # Will cd to this directory, so just pass the filename
        # so the app is not confused by relative paths
        fasta_path = fasta_filename

    index_basename = "%s_%s" % (splitext(fasta_filename)[0], getpid())

    db_name = join(output_dir, index_basename)

    # Instantiate the object
    sdb = IndexDB(WorkingDir=output_dir, HALT_EXEC=HALT_EXEC)

    # The parameter --ref STRING must follow the format where
    # STRING = /path/to/ref.fasta,/path/to/ref.idx
    sdb.Parameters['--ref'].on("%s,%s" % (fasta_path, db_name))

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
    db_filepaths = [v.name for k, v in app_result.items()
                    if k not in {'StdErr', 'StdOut'} and hasattr(v, 'name')]

    return db_name, db_filepaths


class Sortmerna(CommandLineApplication):
    """ SortMeRNA generic application controller for OTU picking
    """

    _command = 'sortmerna'
    _command_delimiter = ' '
    _parameters = {
        # Verbose (log to stdout)
        '-v': FlagParameter('-', Name='v', Value=True),

        # Fasta or Fastq input query sequences file
        '--reads': ValuedParameter('--', Name='reads', Delimiter=' ',
                                   IsPath=True, Value=None),

        # Fasta reference file followed by indexed reference
        '--ref': ValuedParameter('--', Name='ref', Delimiter=' ',
                                 IsPath=True, Value=None),

        # File path + base name for all output files
        '--aligned': ValuedParameter('--', Name='aligned', Delimiter=' ',
                                     IsPath=True, Value=None),

        # Output log file with parameters used to launch sortmerna and
        # statistics on final results (the log file takes on
        # the basename given in --aligned and the extension '.log')
        '--log': FlagParameter('--', Name='log', Value=True),

        # Output Fasta or Fastq file of aligned reads (flag)
        '--fastx': FlagParameter('--', Name='fastx', Value=True),

        # Output BLAST alignment file, options include [0,3] where:
        # 0: Blast-like pairwise alignment,
        # 1: Blast tabular format,
        # 2: 1 + extra column for CIGAR string,
        # 3: 2 + extra column for query coverage
        '--blast': ValuedParameter('--', Name='blast', Delimiter=' ',
                                   IsPath=False, Value=None),

        # Output SAM alignment file
        '--sam': FlagParameter('--', Name='sam', Value=False),

        # Output SQ tags in the SAM file (useful for whole-genome alignment)
        '--SQ': FlagParameter('--', Name='SQ', Value=False),

        # Report the best INT number of alignments
        '--best': ValuedParameter('--', Name='best', Delimiter=' ',
                                  IsPath=False, Value="1"),

        # Report first INT number of alignments
        '--num_alignments': ValuedParameter('--', Name='num_alignments',
                                            Delimiter=' ', IsPath=False,
                                            Value=None),

        # Number of threads
        '-a': ValuedParameter('-', Name='a', Delimiter=' ',
                              IsPath=False, Value="1"),

        # E-value threshold
        '-e': ValuedParameter('-', Name='e', Delimiter=' ',
                              IsPath=False, Value="1"),

        # Similarity threshold
        '--id': ValuedParameter('--', Name='id', Delimiter=' ',
                                IsPath=False, Value="0.97"),

        # Query coverage threshold
        '--coverage': ValuedParameter('--', Name='coverage', Delimiter=' ',
                                      IsPath=False, Value="0.97"),

        # Output Fasta/Fastq file with reads failing to pass the --id and
        # --coverage thresholds for de novo clustering
        '--de_novo_otu': FlagParameter('--', Name='de_novo_otu', Value=True),

        # Output an OTU map
        '--otu_map': FlagParameter('--', Name='otu_map', Value=True),

        # Print a NULL alignment string for non-aligned reads
        '--print_all_reads': FlagParameter('--', Name='print_all_reads',
                                           Value=False)
    }
    _synonyms = {}
    _input_handler = '_input_as_string'
    _supress_stdout = False
    _supress_stderr = False

    def _get_result_paths(self, data):
        """ Set the result paths """

        result = {}

        # get the file extension of the reads file (sortmerna
        # internally outputs all results with this extension)
        fileExtension = splitext(self.Parameters['--reads'].Value)[1]

        # at this point the parameter --aligned should be set as
        # sortmerna will not run without it
        if self.Parameters['--aligned'].isOff():
            raise ValueError("Error: the --aligned parameter must be set.")

        # file base name for aligned reads
        output_base = self.Parameters['--aligned'].Value

        # Blast alignments
        result['BlastAlignments'] =\
            ResultPath(Path=output_base + '.blast',
                       IsWritten=self.Parameters['--blast'].isOn())

        # SAM alignments
        result['SAMAlignments'] =\
            ResultPath(Path=output_base + '.sam',
                       IsWritten=self.Parameters['--sam'].isOn())

        # OTU map (mandatory output)
        result['OtuMap'] =\
            ResultPath(Path=output_base + '_otus.txt',
                       IsWritten=self.Parameters['--otu_map'].isOn())

        # FASTA file of sequences in the OTU map (madatory output)
        result['FastaMatches'] =\
            ResultPath(Path=output_base + fileExtension,
                       IsWritten=self.Parameters['--fastx'].isOn())

        # FASTA file of sequences not in the OTU map (mandatory output)
        result['FastaForDenovo'] =\
            ResultPath(Path=output_base + '_denovo' +
                       fileExtension,
                       IsWritten=self.Parameters['--de_novo_otu'].isOn())
        # Log file
        result['LogFile'] =\
            ResultPath(Path=output_base + '.log',
                       IsWritten=self.Parameters['--log'].isOn())

        return result

    def getHelp(self):
        """Method that points to documentation"""
        help_str = ("SortMeRNA is hosted at:\n"
                    "http://bioinfo.lifl.fr/RNA/sortmerna/\n"
                    "https://github.com/biocore/sortmerna\n\n"
                    "The following paper should be cited if this resource is "
                    "used:\n\n"
                    "Kopylova, E., Noe L. and Touzet, H.,\n"
                    "SortMeRNA: fast and accurate filtering of ribosomal RNAs "
                    "in\n"
                    "metatranscriptomic data, Bioinformatics (2012) 28(24)\n"
                    )
        return help_str


def sortmerna_ref_cluster(seq_path=None,
                          sortmerna_db=None,
                          refseqs_fp=None,
                          result_path=None,
                          tabular=False,
                          max_e_value=1,
                          similarity=0.97,
                          coverage=0.97,
                          threads=1,
                          best=1,
                          aligned='sortmerna_results',
                          HALT_EXEC=False
                          ):
    """Launch sortmerna OTU picker

        Parameters
        ----------
        seq_path : str
            filepath to query sequences.
        sortmerna_db : str
            indexed reference database.
        refseqs_fp : str
            filepath of reference sequences.
        result_path : str
            filepath to output OTU map.
        max_e_value : float, optional
            E-value threshold [default: 1].
        similarity : float, optional
            similarity %id threshold [default: 0.97].
        coverage : float, optional
            query coverage % threshold [default: 0.97].
        threads : int, optional
            number of threads to use (OpenMP) [default: 1].
        tabular : bool, optional
            output BLAST tabular alignments [default: False].
        best : int, optional
            number of best alignments to output per read
            [default: 1]
        aligned : str, optional
            base name for --aligned output files

        Returns
        -------
        clusters : dict of lists
            OTU ids and reads mapping to them

        failures : list
            reads which did not align
    """

    # Instantiate the object
    smr = Sortmerna(HALT_EXEC=HALT_EXEC)

    # Set input query sequences path
    if seq_path is not None:
        smr.Parameters['--reads'].on(seq_path)
    else:
        raise ValueError("Error: a read file is mandatory input.")

    # Set the input reference sequence + indexed database path
    if sortmerna_db is not None:
        smr.Parameters['--ref'].on("%s,%s" % (refseqs_fp, sortmerna_db))
    else:
        raise ValueError("Error: an indexed database for reference set %s must"
                         " already exist.\nUse indexdb_rna to index the"
                         " database." % refseqs_fp)

    if result_path is None:
        raise ValueError("Error: the result path must be set.")

    # Set output results path (for Blast alignments, clusters and failures)
    output_dir = dirname(result_path)
    if output_dir is not None:
        output_file = join(output_dir, aligned)
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

    # Set Blast tabular output
    # The option --blast 3 represents an
    # m8 blast tabular output + two extra
    # columns containing the CIGAR string
    # and the query coverage
    if tabular:
        smr.Parameters['--blast'].on("3")

    # Set number of threads
    if threads is not None:
        smr.Parameters['-a'].on(threads)

    # Run sortmerna
    app_result = smr()

    # Put clusters into a map of lists
    f_otumap = app_result['OtuMap']
    rows = (line.strip().split('\t') for line in f_otumap)
    clusters = {r[0]: r[1:] for r in rows}

    # Put failures into a list
    f_failure = app_result['FastaForDenovo']
    failures = [re.split('>| ', label)[0]
                for label, seq in parse_fasta(f_failure)]

    # remove the aligned FASTA file and failures FASTA file
    # (currently these are re-constructed using pick_rep_set.py
    # further in the OTU-picking pipeline)
    smr_files_to_remove = [app_result['FastaForDenovo'].name,
                           app_result['FastaMatches'].name,
                           app_result['OtuMap'].name]

    return clusters, failures, smr_files_to_remove


def sortmerna_map(seq_path,
                  output_dir,
                  refseqs_fp,
                  sortmerna_db,
                  e_value=1,
                  threads=1,
                  best=None,
                  num_alignments=None,
                  HALT_EXEC=False,
                  output_sam=False,
                  sam_SQ_tags=False,
                  blast_format=3,
                  print_all_reads=True,
                  ):
    """Launch sortmerna mapper

        Parameters
        ----------
        seq_path : str
            filepath to reads.
        output_dir : str
            dirpath to sortmerna output.
        refseqs_fp : str
            filepath of reference sequences.
        sortmerna_db : str
            indexed reference database.
        e_value : float, optional
            E-value threshold [default: 1].
        threads : int, optional
            number of threads to use (OpenMP) [default: 1].
        best : int, optional
            number of best alignments to output per read
            [default: None].
        num_alignments : int, optional
            number of first alignments passing E-value threshold to
            output per read [default: None].
        HALT_EXEC : bool, debugging parameter
            If passed, will exit just before the sortmerna command
            is issued and will print out the command that would
            have been called to stdout [default: False].
        output_sam : bool, optional
            flag to set SAM output format [default: False].
        sam_SQ_tags : bool, optional
            add SQ field to SAM output (if output_SAM is True)
            [default: False].
        blast_format : int, optional
            Output Blast m8 tabular + 2 extra columns for CIGAR
            string and query coverge [default: 3].
        print_all_reads : bool, optional
            output NULL alignments for non-aligned reads
            [default: True].

        Returns
        -------
        dict of result paths set in _get_result_paths()
    """

    if not (blast_format or output_sam):
        raise ValueError("Either Blast or SAM output alignment "
                         "format must be chosen.")

    if (best and num_alignments):
        raise ValueError("Only one of --best or --num_alignments "
                         "options must be chosen.")

    # Instantiate the object
    smr = Sortmerna(HALT_EXEC=HALT_EXEC)

    # Set the input reference sequence + indexed database path
    smr.Parameters['--ref'].on("%s,%s" % (refseqs_fp, sortmerna_db))

    # Set input query sequences path
    smr.Parameters['--reads'].on(seq_path)

    # Set Blast tabular output
    # The option --blast 3 represents an
    # m8 blast tabular output + two extra
    # columns containing the CIGAR string
    # and the query coverage
    if blast_format:
        smr.Parameters['--blast'].on(blast_format)

    # Output alignments in SAM format
    if output_sam:
        smr.Parameters['--sam'].on()
        if sam_SQ_tags:
            smr.Parameters['--SQ'].on()

    # Turn on NULL string alignment output
    if print_all_reads:
        smr.Parameters['--print_all_reads'].on()

    # Set output results path (for Blast alignments and log file)
    output_file = join(output_dir, "sortmerna_map")
    smr.Parameters['--aligned'].on(output_file)

    # Set E-value threshold
    if e_value is not None:
        smr.Parameters['-e'].on(e_value)

    # Set number of best alignments to output per read
    if best is not None:
        smr.Parameters['--best'].on(best)

    # Set number of first alignments passing E-value threshold
    # to output per read
    if num_alignments is not None:
        smr.Parameters['--num_alignments'].on(num_alignments)

    # Set number of threads
    if threads is not None:
        smr.Parameters['-a'].on(threads)

    # Turn off parameters related to OTU-picking
    smr.Parameters['--fastx'].off()
    smr.Parameters['--otu_map'].off()
    smr.Parameters['--de_novo_otu'].off()
    smr.Parameters['--id'].off()
    smr.Parameters['--coverage'].off()

    # Run sortmerna
    app_result = smr()

    return app_result
