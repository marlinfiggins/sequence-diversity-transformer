import Bio.AlignIO.Interfaces
import numpy as np
import pandas as pd
import joblib
import Bio
from Bio import SeqIO, pairwise2
from . import get_metadata


class Parser(object):
    """ A frontend handler and dispatcher for label creation"""

    def __init__(self):
        self.metadata = None
        self.refs = None
        self.data = None
        self.built_data = False  # flag
        pass

    def build_metadata_from_df(self, metadata):
        """
        Get metadata as pd.DataFrame and built as dictionary with id as label.
        :param metadata:
        :return:
        """
        metadata = metadata.set_index("accession_id").to_dict()

        if self.metadata is None:
            self.metadata = metadata
        else:
            assert type(self.metadata) is dict
            self.metadata.update(metadata)

    def build_metadata_from_file(self, metadata_file):
        """
        Get a filename and add the metadata in the frame
        :param metadata_file:
        :return:
        """
        self.build_metadata_from_df(
            get_metadata.load_metadata_from_file(metadata_file)
        )

    def parse_ref_seqs(self):
        """
        Parse the reference sequences into the data structure
        :return:
        """
        #TODO: Change to dictionary to add for easy referencing of sequences
        ref_dict = {}
        for record in Bio.SeqIO.parse(self.ref_data_file, "fasta"):
            ref_dict[record.name] = record # Double check names are correcty parsed in sequence column

        if self.refs is None:
            self.refs = ref_dict
        else:
            assert type(self.refs) is dict
            self.refs.update(ref_dict)

    def build_data(self, data_file):
        """
        Add data safely to avoid key conflicts.
        :param data_file:
        :return:
        """
        tmp_data = Bio.SeqIO.to_dict(SeqIO.parse(data_file, "fasta"))
        # no current data, just add it
        data_dict = {}
        for key in tmp_data.keys():
            new_key = key.split('|')[0]
            data_dict[new_key] = tmp_data[key]

        if self.data is None:
            self.data = data_dict
        else:
            for key in data_dict.keys():
                # Make sure key does not clash
                assert key not in self.data.keys()
                # add key and copy data
                self.data[key] = data_dict[key]

    def parse_data(self):
        """
        Dispatch workers for data and cleanup

        """
        # TODO: Refactor and do MSA on sets with same reference sequences
        with joblib.Parallel(n_jobs=16) as parallel:
            # Execute 1 copy to each of n workers of data
            # Really only need to pass in accession_id's
            worker_results = parallel(joblib.delayed(worker_parser)(
                self.data,
                self.refs,
                self.metadata,
                accession_id) for accession_id in self.metadata.keys)
            res = list(worker_results)  # sync barrier

        acc_list = []
        aln_seqs = []
        scores = []
        labels = []
        for data_tuple in res:
            acc_list.append(data_tuple[0])
            aln_seqs.append(data_tuple[1])
            scores.append(data_tuple[2])
            labels.append(data_tuple[3])
        return {'accession_id': acc_list,
                'aligned_seqs': aln_seqs,
                'scores': scores,
                'labels': labels}


def get_label(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return [ch1 != ch2 for ch1, ch2 in zip(s1, s2)]


def worker_parser(data, refs, metadata, accession_id):
    """
    Compute Hamming distances and return the string rep for closest seq
    and label derived from shortest reference sequence.

    :param data:
    :param refs:
    :param accession_id: string
    :return: tuple
    """

    seq = data[accession_id].seq
    ref_seq = refs[metadata[accession_id]["ref"]].seq
    aligned = pairwise2.align.globalxx(seq, ref_seq)
    aln_seq, aln_ref, score, begin, end = aligned[0]

    return accession_id, aln_seq, aligned.score, get_label(aln_seq, aln_ref)
