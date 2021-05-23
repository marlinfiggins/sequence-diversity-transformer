import numpy as np
import torch
from .tokenizer import tokenize_seq, onehot_tokens, generate_embedding_map
import copy
import math


class Dataset(torch.utils.data.Dataset):
    # embed map should be built once per application and shared between classes/datasets
    def __init__(self, data_triples, token_len: int = 1, dense=False, embed=None, debed=None):
        super().__init__()
        self.full_data = []
        self.data_list = []
        self.label_list = []
        if embed is None and debed is None:
            self.embedding, self.debedding = generate_embedding_map(token_len, dense=dense)
        else:
            self.embedding = embed
            self.debedding = debed

        for data_item in data_triples:
            label_embed = np.asarray(tokenize_seq(data_item[2]))
            data_embed = np.asarray(onehot_tokens(tokenize_seq(data_item[1]), embed_map=self.embedding))
            self.full_data.append({'key': data_item[0],
                                   'seq': data_item[1],
                                   'label_seq': data_item[2],
                                   # Now we tokenize the label and input string to dense numbered vectors
                                   'label': label_embed,
                                   'tokens': data_embed})
            self.data_list.append(data_embed)
            self.label_list.append(label_embed)

        self.start = 0
        self.end = len(self.full_data)

    def embedding_size(self):
        return len(self.embedding.keys())

    def return_dataset_pytorch(self, batch_size=32):
        dataset = torch.utils.data.DataLoader(Dataset=self,  # A neat self-reference!
                                              pin_memory=True, batch_size=batch_size, prefetch_workers=4)
        return dataset

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:  # safegaurd against single-threading
            iter_start = self.start
            iter_end = self.end
        else:  # Split generators across workers
            per_worker = int(math.ceil((self.end - self.start) / float(worker_info.num_workers)))
            worker_id = worker_info.id
            iter_start = self.start + worker_id * per_worker
            iter_end = min(iter_start + per_worker, self.end)
        if self.eval_mode:
            return iter(zip(self.data_list[iter_start:iter_end], self.label_list[iter_start:iter_end]))
        else:
            return iter(self.full_data[iter_start:iter_end])

    def share_embed_debed(self):
        # Returns a deep copy of the embedding and debedding maps
        return copy.deepcopy(self.embedding), copy.deepcopy(self.debedding)

