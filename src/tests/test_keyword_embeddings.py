from src.Core.tokenizer import _generate_keywords_by_len, tokenize_seq, generate_embedding_map, onehot_tokens

def test_gen_keywords_by_len():
    keywords_1 = _generate_keywords_by_len(1)
    assert len(keywords_1) == 4, 'Expected 4 bases on length 1, ACTG'

    keywords_2 = _generate_keywords_by_len(2)
    assert len(keywords_2) == 16, 'Expected 16 bases on length 2, 4 values in position 1, 4 in position 2'

    keywords_3 = _generate_keywords_by_len(3)
    assert len(keywords_2) == 16, 'Expected 16 bases on length 2, 4 values in position 1, 4 in position 2'

