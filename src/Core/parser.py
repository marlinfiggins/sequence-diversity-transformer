import numpy as np
import pandas as pd
import joblib
import Bio
from Bio import SeqIO
from . import get_metadata


class Parser(object):
    """ A frontend handler and dispatcher for label creation"""

    def __init__(self, ref_data_file, metadata_file):
        self.metadata = get_metadata.load_metadata_from_file(metadata_file)
        self.ref_data_file = ref_data_file
        self.ref_labels = []
        self.ref_seqs = []
        self.parse_ref_seqs()
        self.data = None
        self.data_raw = []
        self.labels = []
        self.built_data = False  # flag
        pass

    def parse_ref_seqs(self):
        """
        Parse the reference sequences into the data structure
        :return:
        """
        for record in Bio.SeqIO.parse(self.ref_data_file, "fasta"):
            self.ref_seqs.append(record.seq)
            self.ref_labels.append(record.name)

    def build_data(self, data_file):
        """
        Add data safely to avoid key conflicts.
        :param data_file:
        :return:
        """
        data = Bio.SeqIO.to_dict(SeqIO.parse(data_file, "fasta"))
        # no current data, just add it
        if self.data is None:
            self.data = data
        else:
            for key in data.keys():
                # Make sure key does not clash
                assert key not in self.data.keys()
                # add key and copy data
                self.data[key] = data[key]

    def parse_data(self):
        """ Dispatch workers for data and cleanup"""
        with joblib.Parallel(n_jobs=16) as parallel:
            # Execute 1 copy to each of n workers of data
            worker_results = parallel(joblib.delayed(worker_parser)(self.data, self.ref_seqs, accession_id) for
                                      accession_id in self.metadata.keys())
            res = list(worker_results)  # sync barrier
        acc_list = []
        neighbors = []
        distances = []
        labels = []
        for data_tuple in res:
            acc_list.append(data_tuple[0])
            neighbors.append(data_tuple[1])
            distances.append(data_tuple[2])
            labels.append(data_tuple[3])
        return acc_list, neighbors, distances, labels


def worker_parser(data_dict, refs, accession_id):
    """
    Compute Hamming distances and return the string rep for closest seq
    and label derived from shortest  reference sequence.

    :param data_dict:
    :param refs:
    :param accession_id: string
    :return: tuple
    """
    seq = data_dict[accession_id]
    for ref in refs:
        pass
    raise NotImplemented
    # return accession_id, neighbor, distance, label
