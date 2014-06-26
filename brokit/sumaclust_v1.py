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

from os.path import split, isdir, dirname, isfile, exists, realpath

from skbio.app.util import CommandLineApplication, ResultPath
from skbio.app.parameters import ValuedParameter, FlagParameter


class Sumaclust(CommandLineApplication):
    """ SumaClust generic application controller for de novo OTU picking
    """

    _command = 'sumaclust'
    _command_delimiter = ' '
    _parameters = {
        # Reference sequence length is the shortest
        '-l': FlagParameter('-', Name='l', Value=True),

        # Filepath of the OTU-map
        '-O': ValuedParameter('-', Name='O', Delimiter=' ',
                              Value=None, IsPath=True),

        # Flag '-f' must be passed to deactivate FASTA output
        '-f': FlagParameter('-', Name='f', Value=True),

        # Number of threads
        '-p': ValuedParameter('-', Name='p', Delimiter=' ',
                              Value=1, IsPath=False),

        # Assign sequence to the best matching cluster seed, rather
        # than the first matching cluster (having >= similarity threshold)
        '-e': FlagParameter('-', Name='e', Value=False),

        # Similarity threshold
        '-t': ValuedParameter('-', Name='t', Delimiter=' ',
                              Value=0.97, IsPath=False),

        # Maximum ratio between abundance of two sequences so that the
        # less abundant one can be considered as a variant of the more
        # abundant one.
        '-R': ValuedParameter('-', Name='R', Delimiter=' ',
                              Value=1, IsPath=False)
    }

    _synonyms = {}
    _input_handler = '_input_as_string'
    _supress_stdout = False
    _supress_stderr = False

    def _get_result_paths(self, data):
        """ Set the result paths
        """

        result = {}

        # OTU map (mandatory output)
        result['OtuMap'] = ResultPath(Path=self.Parameters['-O'].Value,
                                      IsWritten=True)

        # SumaClust will not produce any output file if the
        # input file was empty, so we create an empty
        # output file
        if not isfile(result['OtuMap'].Path):
            otumap_f = open(result['OtuMap'].Path, 'w')
            otumap_f.close()

        return result

    def getHelp(self):
        """ Method that points to documentation
        """
        help_str = ("SumaClust is hosted at:\n"
                    "http://metabarcoding.org/sumatra/\n\n"
                    "The following paper should be cited if this resource "
                    "is used:\n\n"
                    "SUMATRA and SUMACLUST: fast and exact comparison and "
                    "clustering "
                    "of full-length barcode sequences\n"
                    "Mercier, C., Boyer, F., Kopylova, E., Taberlet, P., "
                    "Bonin, A. and Coissac E.,"
                    "2014 (in preparation)\n"
                    )

        return help_str


def sumaclust_denovo_cluster(seq_path=None,
                             result_path=None,
                             shortest_len=True,
                             similarity=0.97,
                             threads=1,
                             exact=False,
                             HALT_EXEC=False
                             ):
    """ Function  : launch SumaClust de novo OTU picker

        Parameters: seq_path, filepath to reads;
                    result_path, filepath to output OTU map;
                    shortest_len, boolean;
                    similarity, the similarity threshold (between (0,1]);
                    threads, number of threads to use;
                    exact, boolean to perform exact matching

        Return    : a dictionary of all output files set in
                    _get_result_paths() with the file
                    descriptors pointing to an open file
                    as the values
    """

    # Sequence path is mandatory
    if (seq_path is None
            or not exists(seq_path)):
        raise ValueError("Error: FASTA query sequence filepath is "
                         "mandatory input.")

    # Output directory is mandatory
    if (result_path is None
            or not isdir(dirname(realpath(result_path)))):
        raise ValueError("Error: output directory is mandatory input.")

    # Instantiate the object
    sumaclust = Sumaclust(HALT_EXEC=HALT_EXEC)

    # Set the OTU-map filepath
    sumaclust.Parameters['-O'].on(result_path)

    # Set the similarity threshold
    if similarity is not None:
        sumaclust.Parameters['-t'].on(similarity)

    # Set the option to perform exact clustering (default: False)
    if exact:
        sumaclust.Parameters['-e'].on()

    # Turn off option for reference sequence length to be the shortest
    if not shortest_len:
        sumaclust.Parameters['-l'].off()

    # Set the number of threads
    if threads > 1:
        sumaclust.Parameters['-p'].on(threads)
    elif threads < 1:
        raise ValueError("Number of threads must be positive.")

    # Launch SumaClust,
    # set the data string to include the read filepath
    # (to be passed as final arguments in the sumaclust command)
    app_result = sumaclust(seq_path)

    # Put clusters into a list of lists
    clusters = []
    f_otumap = app_result['OtuMap']
    clusters = [line.strip().split('\t')[1:] for line in f_otumap]

    # Return clusters
    return clusters
