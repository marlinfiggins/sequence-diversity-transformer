from Bio import SeqIO
import pandas as pd


def get_sequence_metadata(file, save=True):
    dicts = {}
    for index, record in enumerate(SeqIO.parse(file, "fasta")):
        split_header = record.description.split("|")

        dicts[index] = {"accession_id": split_header[0],
                        "strain": split_header[1],
                        "year": split_header[2],
                        "month": split_header[3],
                        "day": split_header[4],
                        "segname": split_header[5],
                        "segnum": split_header[6],
                        "length": len(record.seq)}

    dataframe = pd.DataFrame.from_dict(dicts, "index")

    if save:
        dataframe.to_csv(file.split(".fa")[0] + "_metadata.tsv", sep="\t")
    return dataframe


def load_metadata_from_file(filename):
    """
    Helper function to load in metadata from a file.
    :param filename: string
    :return:
    """
    return pd.read_csv(filename, delimiter='\t')


if __name__ == "__main__":
    df = get_sequence_metadata("../data/FASTA_samples_no_identical_seqs.fa" )
