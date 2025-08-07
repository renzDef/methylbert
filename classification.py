''' Second part of the MethylBERT procedure. Fine Tuning needs to be finished. '''

from methylbert.data import finetune_data_generate as fdg
from methylbert.utils import set_seed
from torch.utils.data import DataLoader
from methylbert.data.vocab import MethylVocab
from methylbert.data.dataset import MethylBertFinetuneDataset
from methylbert.trainer import MethylBertFinetuneTrainer
import os
import pandas as pd
from methylbert.deconvolute import deconvolute



''' This part preprocesses the input bulk data used for testing later.
    Here, the output file data.csv is generated.
    This should also be loaded into a dataloader later.

'''

set_seed(42)
seq_len=150
n_mers=3
batch_size=5
num_workers=5

f_bam = "../data/bam_for_classification/sorted_bulk_data.bam"

fdg.finetune_data_generate(
    input_file = f_bam,
    f_dmr = f_dmr,
    f_ref = f_ref,
    output_dir=out_dir,
    n_mers=3, # 3-mer DNA sequences
    n_cores=20
)


# This part is for the actual deconvolution of the bulk data.

output_path="tmp/deconvolution/"


# loading the bulk data
dataset = MethylBertFinetuneDataset("tmp/data.csv",
                                    tokenizer,
                                    seq_len=seq_len)
data_loader = DataLoader(dataset, batch_size=batch_size, num_workers=num_workers)

#restore_dir = "tmp/fine_tune/"
trainer = MethylBertFinetuneTrainer(len(tokenizer),
                                    train_dataloader=data_loader,
                                    test_dataloader=data_loader,
                                    )
trainer.load("model") # Load the fine-tuned MethylBERT model

''' deconvolution produces the following output files:


    decconvolution.csv : tumour deconvolution result
    FI.csv : the Fisher information value
    res.csv : read classification result (the classification result for each read is given in pred column)

    If estimate adjustment is needed, add adjustment = True as parameter to the function.

'''
deconvolute(trainer = trainer,
            tokenizer = tokenizer,
            data_loader = data_loader,
            output_path = output_path,
            df_train = pd.read_csv("tmp/train_seq.csv", sep="\t"))
