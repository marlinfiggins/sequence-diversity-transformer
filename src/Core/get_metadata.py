from Bio import SeqIO
import pandas as pd
import numpy as np


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


# Query function
def query_metadata(df: pd.DataFrame, start_year: int, end_year: int = 2020, strain: str = None,
                   seq_name: str = None, strain_except=None):
    """
    Helper function to filter metadata dataframe.
    Returns a dataframe with rows matching seq_name, strain
    :param strain_except:
    :param strain:
    :param seq_name:
    :param end_year:
    :param start_year:
    :param df: pandas dataframe
           start_year: int
           end_year: (optional) int
           strain: (optional, case sensitive) string
           seq_name: (optional, case sensitive) string
    """
    # Filter by seqname
    if seq_name is not None:
        assert type(seq_name) is not None
        df = df.loc[df['segname'].str.contains(seq_name)]
    # Filter by strain
    if strain is not None:
        assert type(seq_name) is not None
        df = df.loc[df['strain'].str.contains(strain)]
    if strain_except is not None:
        df = df.loc[~df['strain'].str.contains(strain_except)]
    # Filter by year
    years = [str(x) for x in np.arange(start_year, end_year + 1, 1)]
    out = df.loc[df['year'].isin(years)]

    return out


def load_metadata_from_file(filename):
    """
    Helper function to load in metadata from a file.
    :param filename: string
    :return:
    """
    return pd.read_csv(filename, delimiter='\t')
