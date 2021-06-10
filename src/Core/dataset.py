import numpy as np
import torch
from .tokenizer import tokenize_seq, onehot_tokens, generate_embedding_map
import copy
import math


class Dataset(torch.utils.data.IterableDataset):
    # embed map should be built once per application and shared between classes/datasets
    def __init__(self, data_dict, embed): #  token_len: int = 1, dense=False, embed=None, debed=None):
        super().__init__()
        self.full_data = []
        self.data_list = []
        self.label_list = []
        self.embedding = embed
        # if embed is None and debed is None:
        #    self.embedding, self.debedding = generate_embedding_map(token_len, dense=dense)
        # else:
        #    self.embedding = embed
        #    self.debedding = debed

        for idx in range(len(data_dict['aligned_seqs'])):
            # label_embed = np.asarray(tokenize_seq(data_item[2]))
            # data_embed = np.asarray(onehot_tokens(tokenize_seq(data_item[1]), embed_map=self.embedding))

            self.full_data.append({'key': data_dict['accession_id'][idx],
                                   'seq': data_dict['aligned_seqs'][idx],
                                   'label_seq': data_dict['labels'][idx],
                                   # Now we tokenize the label and input string to dense numbered vectors
                                   'tokens': data_dict['numerical_tokens'][idx],
                                   'drop_pos': 0})
            self.data_list.append(torch.tensor(data_dict['numerical_tokens'][idx]))
            self.label_list.append(torch.tensor(data_dict['labels'][idx]))

        self.start = 0
        self.end = len(self.full_data)
        self.eval_mode = True

        self.subseq_split = False
        self.subseq_data = []
        self.subseq_labels = []

    def drop_dashes(self, dash_pos=4, drop_token = '-'):
        for idx, data_item in enumerate(self.full_data):
            # Find the tokens to drop
            drop_pos = np.where(np.asarray(data_item['tokens']) == dash_pos)
            # Store this info
            data_item['drop_pos'] = drop_pos
            # Build a mask
            mask = np.ones(len(data_item['tokens']), dtype=bool)
            mask[drop_pos] = False
            # Drop the tokens
            data_item['tokens'] = np.asarray(data_item['tokens'] )[mask].tolist()
            data_item['label_seq'] = np.asarray(data_item['label_seq'])[mask].tolist()
            data_item['seq'] = [v for v in data_item['seq'] if v != drop_token]
            self.label_list[idx] = torch.as_tensor(self.label_list[idx].numpy()[mask])

    def embedding_size(self):
        return len(self.embedding.keys())

    def make_subseqs(self, length=500):
        for data_entry in self.full_data:
            for idx in range(0, len(data_entry['seq'])-length, 100):
                self.subseq_data.append(torch.as_tensor(data_entry['tokens'][idx:idx+length]))
                self.subseq_labels.append(torch.as_tensor(data_entry['label_seq'][idx:idx+length]))

        self.subseq_split = True
        self.end = len(self.subseq_data)

    def return_dataset_pytorch(self, batch_size=32):
        dataset = torch.utils.data.DataLoader(dataset=self,  # A neat self-reference!
                                              pin_memory=True, batch_size=batch_size)
        return dataset

    def __len__(self):
        return self.end

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
        # Switch which data iterable to look inside
        if self.subseq_split:
            list_data = self.subseq_data
            list_labels = self.subseq_labels
        else:
            list_data = self.full_data
            list_labels = self.label_list

        if self.eval_mode:
            return iter(zip(list_data[iter_start:iter_end], list_labels[iter_start:iter_end]))
        else:
            return iter(self.data_list[iter_start:iter_end])

    def share_embed(self):
        # Returns a deep copy of the embedding map
        return copy.deepcopy(self.embedding)

