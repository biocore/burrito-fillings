#!/opt/python-2.7.3/bin python

#-----------------------------------------------------------------------------
# Copyright (c) 2015--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

""" Application controller for vsearch v1.0.7 """

from os.path import splitext, abspath, join, dirname
from tempfile import mkstemp
from shutil import rmtree

from skbio.parse.sequences import parse_fasta
from burrito.parameters import ValuedParameter, FlagParameter
from burrito.util import (CommandLineApplication, ResultPath,
                          ApplicationError, ApplicationNotFoundError)
from skbio.util import remove_files


class Vsearch(CommandLineApplication):

    """ Vsearch ApplicationController """

    _command = 'vsearch'
    _input_handler = '_input_as_parameters'
    _parameters = {
        # Output to specified FASTA file
        '--output': ValuedParameter('--', Name='output', Delimiter=' ',
                                    IsPath=True),

        # Filename for UCLUST-like output
        '--uc': ValuedParameter('--', Name='uc', Delimiter=' ',
                                IsPath=True),
        
        # Filename for BLAST-like tab-separated output
        '--blast6out': ValuedParameter('--', Name='blast6out', Delimiter=' ',
                                       IsPath=True),

        # ID percent for OTU, by default is 97%
        '--id': ValuedParameter('--', Name='id', Delimiter=' ',
                                IsPath=False, Value="0.97"),
        
        # ID definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)
        '--iddef': ValuedParameter('--', Name='iddef',
                                   Delimiter=' ', IsPath=False,
                                   Value="2"),

        # Number of hits to accept and show per strand (1)
        '--maxaccepts':
        ValuedParameter('--', Name='maxaccepts', Delimiter=' ', Value="1"),

        # Number of non-matching hits to consider (32)
        '--maxrejects':
        ValuedParameter('--', Name='maxrejects', Delimiter=' ', Value="32"),

        # Indicate that input sequences are presorted
        '--usersort': FlagParameter('--', Name='usersort'),

        # Take into account the abundance annotations present
        # in the input fasta file
        '--sizein' : FlagParameter('--', Name='sizein'),

        # Add abundance annotations to the output fasta files
        '--sizeout': FlagParameter('--', Name='sizeout'),

        # Dereplicate exact sequences in the given FASTA file
        '--derep_fulllength': ValuedParameter('--', Name='derep_fulllength',
                                              Delimiter=' ', IsPath=True),

        # Dereplicate plus or both strands (both)
        '--strand': ValuedParameter('--', Name='strand', Delimiter=' ',
                                    IsPath=False, Value="both")

        # Discard sequences with an abundance value greater than integer
        '--maxuniquesize': ValuedParameter('--', Name='maxuniquesize', Delimiter=' ',
                                           IsPath=False),

        # Discard sequences with an abundance value smaller than integer
        '--minuniquesize': ValuedParameter('--', Name='minuniquesize', Delimiter=' ',
                                           IsPath=False),

        # Abundance sort sequences in given FASTA file
        '--sortbysize': ValuedParameter('--', Name='sortbysize', Delimiter=' ',
                                      IsPath=True),
        
        # When using --sortbysize, discard sequences
        # with an abundance value greater than maxsize
        '--maxsize': ValuedParameter('--', Name='maxsize', Delimiter=' ', IsPath=False),

        # When using --sortbysize, discard sequences
        # with an abundance value smaller than minsize
        '--misize': ValuedParameter('--', Name='minsize', Delimiter=' ', IsPath=False),

        # Output cluster consensus sequences to FASTA file
        '--consout': ValuedParameter('--', Name='consout', Delimiter=' ',
                                     IsPath=True),

        # Chimera detection: min abundance ratio of parent vs chimera (2.0)
        '--abskew': ValuedParameter('--', Name='abskew', Delimiter=' ',
                                    IsPath=False, Value="2.0"),
        # Detect chimeras de novo
        '--uchime_denovo': ValuedParameter('--', Name='uchime_denovo', Delimiter=' ',
                                    IsPath=True),
        
        # Detect chimeras using a reference database
        '--uchime_ref': ValuedParameter('--', Name='uchime_ref',
                                        Delimiter=' ', IsPath=True)

        # Output chimera alignments to 3-way alignment file (filepath)
        '--uchimealns': ValuedParameter('--', Name='uchimealns', Delimiter=' ',
                                      IsPath=True),

        # Output chimeric sequences to file (filepath)
        '--chimeras': ValuedParameter('--', Name='chimeras',
                                      Delimiter=' ', IsPath=True)

        # Output non-chimera filepath
        '--nonchimeras': ValuedParameter('--', Name='nonchimeras',
                                         Delimiter=' ', IsPath=True),

        # Reference database for --uchime_ref
        '--db': ValuedParameter('--', Name='db', Delimiter=' ', IsPath=True),

        # Output to chimera info to tab-separated file
        '--uchimeout': ValuedParameter('--', Name='uchimeout', Delimiter=' ',
                                       IsPath=True),
        
        # Number of computation threads to use (1 to 256)
        '--threads': ValuedParameter('--', Name='threads', Delimiter=' ',
                                     IsPath=False, Value="1")
    }

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['--uc', '--output', '--sortbysize',
                          '--consout', '--uchime_denovo',
                          '--derep_fulllength', '--maxuniquesize',
                          '--minuniquesize', '--sizein',
                          '--sizeout', '--strand', '--threads',
                          '--uchime_ref', '--chimeras',
                          '--nonchimeras', '--db', '--uchimeout',
                          '--blast6out', '--abskew',
                          '--sortbysize', '--maxsize', '--minsize']

        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError(
                "Unsupported parameter(s) passed when calling vsearch: %s" %
                ' '.join(unsupported_parameters))

        for v in allowed_values:
            # turn the parameter off so subsequent runs are not
            # affected by parameter settings from previous runs
            self.Parameters[v].off()
            if v in data:
                # turn the parameter on if specified by the user
                self.Parameters[v].on(data[v])

        return ''

    def _get_result_paths(self, data):
        """ Set the result paths """

        result = {}

        result['Output'] = ResultPath(
            Path=self.Parameters['--output'].Value,
            IsWritten=self.Parameters['--output'].isOn())

        result['ClusterFile'] = ResultPath(
            Path=self.Parameters['--uc'].Value,
            IsWritten=self.Parameters['--uc'].isOn())

        # uchime 3-way global alignments
        result['Output_aln'] = ResultPath(
            Path=self.Parameters['--uchimealns'].Value,
            IsWritten=self.Parameters['--uchimealns'].isOn())

        # uchime tab-separated format
        result['Output_tabular'] = ResultPath(
            Path=self.Parameters['--uchimeout'].Value,
            IsWritten=self.Parameters['--uchimeout'].isOn())

        # chimeras fasta file output
        result['Output_chimeras'] = ResultPath(
            Path=self.Parameters['--chimeras'].Value,
            IsWritten=self.Parameters['--chimeras'].isOn())

        # nonchimeras fasta file output
        result['Output_nonchimeras'] = ResultPath(
            Path=self.Parameters['--nonchimeras'].Value,
            IsWritten=self.Parameters['--nonchimeras'].isOn())

        return result


    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        VSEARCH is hosted at:
        https://github.com/torognes/vsearch

        The following papers should be cited if this resource is used:

        Paper pending. Please cite the github page in the meanwhile. 
        """
        return help_str


def vsearch_dereplicate_exact_seqs(
        fasta_filepath,
        output_filepath,
        output_uc=False,
        working_dir=None,
        strand="both",
        maxuniquesize=None,
        minuniquesize=None,
        sizein=False,
        sizeout=True,
        log_name="derep.log",
        HALT_EXEC=False):
    """ Generates clusters and fasta file of
        dereplicated subsequences

        Parameters
        ----------
        
        fasta_filepath : string
           input filepath of fasta file to be dereplicated
        output_filepath : string
           write the dereplicated sequences to output_filepath
        working_dir : string, optional
           directory path for storing intermediate output
        output_uc : boolean, optional
           uutput dereplication results in a file using a
           uclust-like format
        strand : string, optional
           when searching for strictly identical sequences,
           check the 'strand' only (default: both) or
           check the plus strand only
        maxuniquesize : integer, optional
           discard sequences with an abundance value greater
           than maxuniquesize
        minuniquesize : integer, optional
           discard sequences with an abundance value smaller
           than integer
        sizein : boolean, optional
           take into account the abundance annotations present in
           the input fasta file,  (search for the pattern
           "[>;]size=integer[;]" in sequence headers)
        sizeout : boolean, optional
           add abundance annotations to the output fasta file
           (add the pattern ";size=integer;" to sequence headers)
        log_name : string, optional
           specifies log filename
        HALT_EXEC : boolean, optional
           used for debugging app controller
    """

    # write all vsearch output files to same directory
    # as output_filepath if working_dir is not specified
    if not working_dir:
        working_dir = dirname(abspath(output_filepath))

    app = Vsearch(WorkingDir=working_dir, HALT_EXEC=HALT_EXEC)

    log_filepath = join(working_dir, log_name)
    uc_filepath = None
    if output_uc:
        uc_filepath = join(working_dir, 'vsearch_uc_dereplicated.uc')
        app.Parameters['--uc'].on(uc_filepath)

    if maxuniquesize:
        app.Parameters['--maxuniquesize'].on(maxuniquesize)
    if minuniquesize:
        app.Parameters['--minuniquesize'].on(minuniquesize)
    if sizein:
        app.Parameters['--sizein'].on(sizein)
    if sizeout:
        app.Parameters['--sizeout'].on(sizeout)
    if (strand == "both" or strand == "plus"):
        app.Parameters['--strand'].on(strand)
    else:
        raise ValueError("Option --strand accepts only 'both'
                          or 'plus' values")
    app.Parameters['--derep_fulllength'].on(fasta_filepath)
    app.Parameters['--output'].on(output_filepath)

    app_result = app()

    return app_result, output_filepath, uc_filepath


def vsearch_sort_by_abundance(
        fasta_filepath,
        output_filepath,
        working_dir=None,
        minsize=None,
        maxsize=None,
        log_name="abundance_sort.log",
        HALT_EXEC=False):
    """ Fasta entries are sorted by decreasing abundance
        (Fasta entries are assumed to be dereplicated with
        the pattern "[>;]size=integer[;]" present in the
        read label, ex. use function vsearch_dereplicate_exact_seqs
        prior to calling this function)

        Parameters
        ----------

        fasta_filepath : string
           input fasta file (dereplicated fasta)
        output_filepath : string
           output filepath for the sorted sequences in fasta format
        working_dir : string, optional
           working directory to store intermediate files
        minsize : integer, optional
           discard sequences with an abundance value smaller than
           minsize
        maxsize : integer, optional
           discard sequences with an abundance value greater than
           maxsize
        log_name : string, optional
           log filename
        HALT_EXEC : boolean, optional
           used for debugging app controller
    """

    # set working dir to same directory as the output
    # file (if not set)
    if not working_dir:
        working_dir = os.path.dirname(output_filepath) 

    app = Vsearch(WorkingDir=working_dir, HALT_EXEC=HALT_EXEC)

    log_filepath = join(working_dir, log_name)

    if minsize:
        app.Parameters['--minsize'].on(minsize)

    if maxsize:
        app.Parameters['--maxsize'].on(maxsize)

    app.Parameters['--sortbysize'].on(fasta_filepath)
    app.Parameters['--output'].on(output_filepath)

    app_result = app()

    return app_result, output_filepath


def vsearch_chimera_filter_de_novo(
        fasta_filepath,
        working_dir,
        output_chimeras=True,
        output_nonchimeras=True,
        output_alns=False,
        output_tabular=False,
        log_name="vsearch_uchime_de_novo_chimera_filtering.log",
        HALT_EXEC=False):
    """ Detect chimeras present in the fasta-formatted filename,
        without external references (i.e. de novo). Automatically
        sort the sequences in filename by decreasing abundance
        beforehand. Output chimeras and non-chimeras to FASTA files
        and/or 3-way global alignments and/or tabular output.

        Parameters
        ----------

        fasta_filepath : string
           input fasta file (dereplicated fasta)
        working_dir : string
           directory path for all output files
        output_chimeras : boolean, optional
           output chimeric sequences to file, in fasta format
        output_nonchimeras : boolean, optional
           output nonchimeric sequences to file, in fasta format
        output_alns : boolean, optional
           output 3-way global alignments (parentA, parentB, chimera)
           in human readable format to file
        output_tabular : boolean, optional
           output results using the uchime tab-separated format of
           18 fields (see Vsearch user manual)
        HALT_EXEC : boolean, optional
           used for debugging app controller
    """

    app = Vsearch(WorkingDir=working_dir, HALT_EXEC=HALT_EXEC)

    if not (output_chimeras or
            output_nonchimeras or
            output_alns or
            output_tabular):
        raise ValueError("At least one output format (output_chimeras,
                          output_nonchimeras, output_alns, output_tabular)
                          must be selected")

    output_chimera_filepath = None
    output_non_chimera_filepath = None
    output_alns_filepath = None
    output_tabular_filepath = None

    # set output filepaths
    if output_chimeras:
        output_chimera_filepath = join(working_dir, 'uchime_chimeras.fasta')
        app.Parameters['--chimeras'].on(output_chimera_filepath)
    if output_nonchimeras:
        output_non_chimera_filepath = join(working_dir, 'uchime_non_chimeras.fasta')
        app.Parameters['--nonchimeras'].on(output_non_chimera_filepath)
    if output_alns:
        output_alns_filepath = join(working_dir, 'uchime_alignments.txt')
        app.Parameters['--uchimealns'].on(output_alns_filepath)
    if output_tabular:
        output_tabular_filepath = join(working_dir, 'uchime_tabular.txt')
        app.Parameters['--uchimeout'].on(output_tabular_filepath)
    
    app.Parameters['--uchime_denovo'].on(fasta_filepath)

    log_filepath = join(working_dir, log_name)

    app_result = app()

    return app_result, output_chimera_filepath,
           output_non_chimera_filepath, output_alns_filepath,
           output_tabular_filepath


def vsearch_chimera_filter_ref(
        fasta_filepath,
        working_dir,
        db_filepath,
        output_chimeras=True,
        output_nonchimeras=True,
        output_alns=False,
        output_tabular=False,
        log_name="vsearch_uchime_ref_chimera_filtering.log",
        threads=1,
        HALT_EXEC=False):
    """ Detect chimeras present in the fasta-formatted filename,
        with an external reference (i.e. database). Output
        chimeras and non-chimeras to FASTA files and/or 3-way
        global alignments and/or tabular output.

        Parameters
        ----------

        fasta_filepath : string
           input fasta file (dereplicated fasta)
        working_dir : string
           directory path for all output files
        db_filepath : string
           filepath to reference database
        output_chimeras : boolean, optional
           output chimeric sequences to file, in fasta format
        output_nonchimeras : boolean, optional
           output nonchimeric sequences to file, in fasta format
        output_alns : boolean, optional
           output 3-way global alignments (parentA, parentB, chimera)
           in human readable format to file
        output_tabular : boolean, optional
           output results using the uchime tab-separated format of
           18 fields (see Vsearch user manual)
        threads : integer, optional
           number of computation threads to use (1 to 256)
        HALT_EXEC : boolean, optional
           used for debugging app controller
    """

    app = Vsearch(WorkingDir=working_dir, HALT_EXEC=HALT_EXEC)

    if not (output_chimeras or
            output_nonchimeras or
            output_alns or
            output_tabular):
        raise ValueError("At least one output format (output_chimeras,
                          output_nonchimeras, output_alns, output_tabular)
                          must be selected")

    output_chimera_filepath = None
    output_non_chimera_filepath = None
    output_alns_filepath = None
    output_tabular_filepath = None

    # set output filepaths
    if output_chimeras:
        output_chimera_filepath = join(working_dir, 'uchime_chimeras.fasta')
        app.Parameters['--chimeras'].on(output_chimera_filepath)
    if output_nonchimeras:
        output_non_chimera_filepath = join(working_dir, 'uchime_non_chimeras.fasta')
        app.Parameters['--nonchimeras'].on(output_non_chimera_filepath)
    if output_alns:
        output_alns_filepath = join(working_dir, 'uchime_alignments.txt')
        app.Parameters['--uchimealns'].on(output_alns_filepath)
    if output_tabular:
        output_tabular_filepath = join(working_dir, 'uchime_tabular.txt')
        app.Parameters['--uchimeout'].on(output_tabular_filepath)

    app.Parameters['--db'].on(db_filepath)
    app.Parameters['--uchime_ref'].on(fasta_filepath)
    
    log_filepath = join(working_dir, log_name)

    app_result = app()

    return app_result, output_chimera_filepath,
           output_non_chimera_filepath, output_alns_filepath,
           output_tabular_filepath

