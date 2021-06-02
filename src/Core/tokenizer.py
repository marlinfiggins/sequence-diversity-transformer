from Bio import Seq
from collections.abc import Callable
from typing import Iterable, Optional


def tokenize_seq(seq: Iterable, token_len: int = 1, pattern: str = 'overlap', len_rule='strict') -> list[str]:
    """
    Tokenize a string by character strings of length `token_len`.
    If parameter token_len == 1`, `pattern` is ignored.
    If `pattern='overlap'` and `parameter token_len > 1`, then characters are repeated between tokens.
    If `pattern='sparse'` characters are not repeated between tokens.
    Parameter `len_rule` controls how to handle the final token. Ignored if `pattern='overlap'` or if
    `token_len = 1`. Default option `'strict'` drops the final codon.
    Option `'weak'` returns an incomplete codon.

    :param len_rule: String - control clipping behavior for final element. Ignored if token_len > 1
    :param seq:
    :param token_len:
    :param pattern:
    :return: list of tokens
    """
    tokens = []
    if pattern == 'overlap' or token_len == 1:
        for idx in range(len(seq) - token_len + 1):
            # Start at the begining
            tokens.append(str(seq[idx:idx + token_len]))
        if token_len > 1:
            for idx in range(token_len-1, 0, -1):
                tokens.append(str(seq[-idx:]))

    elif pattern == 'sparse':
        num_expected_tokens = len(seq) // token_len
        if len_rule == 'strict':
            pass
        elif len_rule == 'weak':
            num_expected_tokens += 1

        for idx in range(num_expected_tokens):
            tokens.append(seq[idx * token_len: (idx + 1) * token_len])
    else:
        raise ValueError(f"Value of {pattern} was not 'overlap' or 'sparse'")
    return tokens


def _generate_keywords_by_len(length: int = 1) -> (list[str], list[str]):
    """

    :param length: int >= 1
    :return: (List of all strings, generated as an intermediate, list of string of length `length`)
    """
    if length > 5:
        raise Warning('Are you sure you want to do that? The combinatorics are not on your side')
    keywords = list('ACTG-')  # str is an iterable!
    keywords_last = []
    exit_flag = False
    for item in keywords:
        for char in 'ACTG-':
            new_item = item + char
            if len(new_item) == length:
                keywords_last.append(new_item)
            if len(new_item) > length:
                exit_flag = True
                break
            keywords.append(new_item)
        if exit_flag:
            break
    return keywords, keywords_last


def generate_embedding_map(seq_len: int = 0, dense=False) -> (Callable[[str], int], Callable[[int], str]):
    """
    Return dicts (mappings) between sequences (str) and dense numberings.
    If dense, uses all sequences of length `<= seq_len`.
    :param dense: Bool,
    :param seq_len:
    :return:
    """

    keywords, keywords_short = _generate_keywords_by_len(seq_len)
    if not dense:  # Only have the last few (full sequences), ignore the shorter sequences
        keywords = keywords_short  # dump short seqs of full values only.

    discrete_embedding = {}  # convert seqs to numbers
    discrete_debedding = {}
    for idx, seq in enumerate(keywords):
        discrete_embedding[seq] = idx
        discrete_debedding[idx] = seq

    return discrete_embedding, discrete_debedding


def onehot_tokens(token_list: list[str], min_token_len: int = 1, dense=True, embed_map: Callable[[str], int] = None) ->\
        (list[int], (Callable[[str], int], Optional[Callable[[int], str]])):
    """
    Top level function to convert a list of tokens into a list of numbers. Returns a triple of
    (numerical_tokens, embedding map, debedding map)
    If no embed function is specified, one is automatically constructed. The associated inverse map is returned as the
    debedding map if the embedding is automatically constructed. If not constructed, `None` is returned instead.
    also returned.
    :param token_list:
    :param min_token_len:
    :param dense:
    :param embed_map:
    :return:
    """

    debed_map = None  # Allocate this variable name now for later possible return
    if embed_map is None:
        # we need to construct the embedding map
        # Figure out how many base pairs per token our map needs to represent
        num_tokens = max(1, min_token_len, max([len(token) for token in token_list]))
        embed_map, debed_map = generate_embedding_map(num_tokens, dense=dense)

    token_numerical = [embed_map[token] for token in token_list]

    return token_numerical, embed_map, debed_map

