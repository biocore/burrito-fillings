"""
Application controller for Swarm version 1.2.7
==============================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, brokit development team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import split, isdir, exists, join
from tempfile import mkstemp
from os import close
from subprocess import Popen, PIPE, STDOUT
import re

from skbio.app.util import CommandLineApplication, ResultPath
from skbio.app.parameters import ValuedParameter, FlagParameter
from skbio.parse.sequences import parse_fasta
from skbio.util.misc import remove_files


class Swarm(CommandLineApplication):
    """ Swarm generic application controller for de novo OTU picking
    """

    _command = 'swarm'
    _command_delimiter = ' '
    _parameters = {
        # Resolution
        '-d': ValuedParameter('-', Name='d', Delimiter=' ',
                              Value=1, IsPath=False),
        # OTU-map result filename
        '-o': ValuedParameter('-', Name='o', Delimiter=' ',
                              Value=None, IsPath=True),
        # Threads
        '-t': ValuedParameter('-', Name='t', Delimiter=' ',
                              Value=1, IsPath=False),
    }

    _synonyms = {}
    _input_handler = '_input_as_string'
    _supress_stdout = False
    _supress_stderr = False
    files_to_remove = []

    def _swarm_breaker(self, seq_path=None):
        """
            Input : a filepath to de-replicated
                    input FASTA reads

            Method: using swarm_breaker.py, break
                    chains of amplicons based on
                    abundance information. Abundance
                    is stored after the final
                    underscore '_' in each sequence
                    label (recommended procedure for
                    Swarm)

            Return: clusters, a list of lists
        """
        swarm_breaker_command = ["swarm_breaker.py",
                                 "-f",
                                 # de-replicated FASTA sequences
                                 seq_path,
                                 "-s",
                                 # Swarm OTU-map
                                 self.Parameters['-o'].Value,
                                 "-d",
                                 # resolution
                                 str(self.Parameters['-d'].Value)]

        try:
            # launch swarm_breaker.py as a subprocess,
            # pipe refined OTU-map to stdout
            proc = Popen(swarm_breaker_command,
                         stderr=STDOUT,
                         stdout=PIPE,
                         close_fds=True)

            # store refined clusters in list of lists
            clusters = []
            for line in proc.stdout:
                seq_ids = re.split("\t| ", line.strip())
                # remove the abundance information from the labels
                for i in range(len(seq_ids)):
                    amplicon, abundance = seq_ids[i].rsplit("_", 1)
                    seq_ids[i] = amplicon
                clusters.append(seq_ids)
        except (OSError, 2):
            raise OSError("Cannot find swarm_breaker.py "
                          "in the $PATH directories.\n")
        proc.wait()

        return clusters

    def _prefilter_exact_matches(self, seqs):
        """
        """
        unique_sequences = {}
        seq_id_map = {}
        filtered_seqs = []
        for seq_id, seq in seqs:
            seq_id = seq_id.split()[0]
            try:
                temp_seq_id = unique_sequences[seq]
            except KeyError:
                temp_seq_id = 'QiimeExactMatch.%s' % seq_id
                unique_sequences[seq] = temp_seq_id
                seq_id_map[temp_seq_id] = []
                filtered_seqs.append((temp_seq_id, seq))
            seq_id_map[temp_seq_id].append(seq_id)
        return filtered_seqs, seq_id_map

    def _apply_identical_sequences_prefilter(self, seq_path):
        """
            Input : a filepath to input FASTA reads
            Method: prepares and writes de-replicated reads
                    to a temporary FASTA file, calls
                    parent method to do the actual
                    de-replication
            Return: exact_match_id_map, a dictionary storing
                    de-replicated amplicon ID as key and
                    all original FASTA IDs with identical
                    sequences as values;
                    unique_seqs_fp, filepath to FASTA file
                    holding only de-replicated sequences
        """
        # creating mapping for de-replicated reads
        seqs_to_cluster, exact_match_id_map =\
            self._prefilter_exact_matches(parse_fasta(open(seq_path, 'U')))

        # create temporary file for storing the de-replicated reads
        fd, unique_seqs_fp = mkstemp(
            prefix='SwarmExactMatchFilter', suffix='.fasta')
        close(fd)

        self.files_to_remove.append(unique_seqs_fp)

        # write de-replicated reads to file
        unique_seqs_f = open(unique_seqs_fp, 'w')
        for seq_id, seq in seqs_to_cluster:
            unique_seqs_f.write('>%s_%d\n%s\n'
                                % (seq_id,
                                   len(exact_match_id_map[seq_id]),
                                   seq))
        unique_seqs_f.close()
        del(seqs_to_cluster)

        return exact_match_id_map, unique_seqs_fp

    def _map_filtered_clusters_to_full_clusters(self, clusters, filter_map):
        """
            Input:  clusters, a list of cluster lists
                    filter_map, the seq_id in each clusters
                                is the key to the filter_map
                                containing all seq_ids with
                                duplicate FASTA sequences
            Output: an extended list of cluster lists
        """
        results = []
        for cluster in clusters:
            full_cluster = []
            for seq_id in cluster:
                full_cluster += filter_map[seq_id]
            results.append(full_cluster)
        return results

    def _get_result_paths(self, data):
        """ Set the result paths
        """

        result = {}

        # Swarm OTU map (mandatory output)
        result['OtuMap'] = ResultPath(Path=self.Parameters['-o'].Value,
                                      IsWritten=True)

        return result

    def getHelp(self):
        """ Method that points to documentation
        """
        help_str = ("Swarm is hosted at:\n"
                    "https://github.com/torognes/swarm\n\n"
                    "The following paper should be cited if this resource "
                    "is used:\n\n"
                    "Swarm: robust and fast clustering method for "
                    "amplicon-based studies\n"
                    "Mahe, F., Rognes, T., Quince, C., de Vargas, C., "
                    "and Dunthorn, M."
                    "2014 (submitted)\n"
                    )

        return help_str


def swarm_denovo_cluster(seq_path=None,
                         output_dir=None,
                         d=1,
                         threads=1,
                         prefilter_identical_sequences=True,
                         HALT_EXEC=False):
    """ Function  : launch the Swarm de novo OTU picker

        Parameters: seq_path, filepath to reads
                    output_dir, output directory path
                    d, resolution
                    threads, number of threads to use
                    prefilter_identical_sequences, boolean
                    whether to dereplicate the sequences prior
                    to clustering them

        Return    : clusters, list of lists
    """

    # Sequence path is mandatory
    if (seq_path is None
            or not exists(seq_path)):
        raise ValueError("Error: FASTA query sequence filepath is "
                         "mandatory input.")

    # Output directory is mandatory
    if (output_dir is None
            or not isdir(output_dir)):
        raise ValueError("Error: output directory is mandatory input.")

    # Instantiate the object
    swarm = Swarm(HALT_EXEC=HALT_EXEC)

    # Set the resolution
    if d > 0:
        swarm.Parameters['-d'].on(d)
    else:
        raise ValueError("Resolution -d must be a positive integer.")

    # Set the number of threads
    if threads > 0:
        swarm.Parameters['-t'].on(threads)
    else:
        raise ValueError("Number of threads must be a positive integer.")

    # Set the result path for Swarm output
    swarm_OTU_map = join(output_dir, "swarm_otus_tmp.txt")
    swarm.Parameters['-o'].on(swarm_OTU_map)

    # Remove this file later, the final OTU-map
    # is output by swarm_breaker.py
    swarm.files_to_remove.append(swarm_OTU_map)

    # De-replicate query sequences
    if prefilter_identical_sequences:
        exact_match_id_map, seq_path =\
            swarm._apply_identical_sequences_prefilter(seq_path)

    # Launch Swarm
    # set the data string to include the read filepath
    # (to be passed as final arguments in the swarm command)
    app_result = swarm(seq_path)

    # Launch swarm_breaker.py to refine Swarm OTUs
    # if the sequences were not de-replicated using QIIME,
    # verify that the labels include abundance after the
    # last occurrence of '_'
    # (this is required for swarm_breaker)
    if not prefilter_identical_sequences:
        with open(seq_path) as file:
            for label, seq in parse_fasta(file):
                try:
                    amplicon, abundance = label.strip("\n").rsplit("_", 1)
                    break
                except ValueError:
                    raise ValueError("FASTA labels must have the abundance "
                                     "information after the last occurrence "
                                     "of separator '_' in the label, eg. "
                                     ">sampleID_sequenceID_abundance")

    # Run swarm_breaker.py to refine the clusters
    clusters = swarm._swarm_breaker(seq_path)

    # Expand clusters
    if prefilter_identical_sequences:
        clusters = swarm._map_filtered_clusters_to_full_clusters(
            clusters, exact_match_id_map)

    remove_files(swarm.files_to_remove, error_on_missing=False)

    # Return clusters
    return clusters
