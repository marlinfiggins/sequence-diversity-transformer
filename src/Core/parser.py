import Bio.AlignIO.Interfaces
import numpy as np
import pandas as pd
import joblib
import Bio
from Bio import SeqIO, pairwise2
from . import get_metadata
import multiprocessing
import tqdm


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
            print('Setting Metadata attribute')
            self.metadata = metadata
        else:
            assert type(self.metadata) is dict
            for key in self.metadata.keys():
                # update each attribute dict
                self.metadata[key].update(metadata[key])

    def build_metadata_from_file(self, metadata_file):
        """
        Get a filename and add the metadata in the frame
        :param metadata_file:
        :return:
        """
        self.build_metadata_from_df(
            get_metadata.load_metadata_from_file(metadata_file)
        )

    def parse_ref_seqs(self, ref_data_file):
        """
        Parse the reference sequences into the data structure
        :return:
        """
        # TODO: Change to dictionary to add for easy referencing of sequences
        ref_dict = {}
        for record in Bio.SeqIO.parse(ref_data_file, "fasta"):
            ref_dict[record.name] = record  # Double check names are correcty parsed in sequence column
            ref_dict[record.name].key_type = '4 (HA)' if 'hemagglutinin' in ref_dict[
                record.name].description.lower() else '6 (NA)'
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

    def parse_data(self, mc=True, save=False, file_name: str = ''):
        """
        Dispatch workers for data and cleanup

        """
        # TODO: Refactor and do MSA on sets with same reference sequences
        with joblib.Parallel(n_jobs=int(multiprocessing.cpu_count() / 2), verbose=10, pre_dispatch='all') as parallel:
            # Execute 1 copy to each of n workers of data
            # Really only need to pass in accession_id's
            if mc:
                worker_results = parallel(joblib.delayed(worker_parser)
                                          ((self.data[accession_id].seq, self.refs, self.metadata, accession_id))
                                          for accession_id in self.metadata['strain'].keys())
            else:
                worker_results = [worker_parser((self.data[accession_id].seq, self.refs, self.metadata, accession_id))
                                  for accession_id in tqdm.tqdm(self.metadata['strain'].keys())]
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
        data_dict = {'accession_id': acc_list,
                     'aligned_seqs': aln_seqs,
                     'scores': scores,
                     'labels': labels}
        return data_dict


def get_label(s1, s2):
    """
    Compute a boolean label between two strings, s1 and s2. Returns a np.ndarray of
    0 and 1 values.
    :param s1: String one
    :param s2: String 2
    :return: np.ndarray of 0's and 1's.
    """
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return np.asarray([ch1 != ch2 for ch1, ch2 in zip(s1, s2)], dtype=int)


def worker_parser(packed_tuple: tuple):
    """
    Compute alignments and return the string rep for closest seq
    and label derived from shortest reference sequence.

    :param packed_tuple: Data, refs, accession_id
    :return: tuple
    """
    seq, refs, metadata, accession_id = packed_tuple
    #seq = data[accession_id].seq
    aligned_best_score = np.inf
    cur_key = 'None'
    al_best = None
    # Todo implement lookup for ref seq
    idx = 0
    for key in refs.keys():
        ref_seq = refs[key].seq
        if refs[key].key_type != metadata['segname'][accession_id]:
            # print(metadata['segname'][accession_id])
            # print(refs[key].key_type)
            # Keys cannot match, skip cases
            # print('skipping case', refs[key].key_type, '!=', metadata['segname'][accession_id])
            continue
        al = pairwise2.align.globalxx(seq, refs[key].seq, one_alignment_only=True)
        if al[0].score < aligned_best_score:
            aligned_best_score = al[0].score
            cur_key = key
            al_best = al

    # ref_seq = refs[metadata[accession_id]["ref"]].seq
    # aligned = pairwise2.align.globalxx(seq, ref_seq)
    aln_seq, aln_ref, score, start, end = al_best[idx].seqA, al_best[idx].seqB, al_best[idx].score, al_best[idx].start, \
                                          al_best[idx].end

    label = get_label(al_best[idx].seqA, al_best[idx].seqB)
    return accession_id, aln_seq, aligned_best_score, label
