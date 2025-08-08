''' Fine tuning of the methylbert model. '''


''' Input data needed:

    Input bulk sample as a BAM/SAM file
    DMRs as a tab-separated .csv file
    Pure tumour and normal samples as BAM/SAM files

'''
from methylbert.data import finetune_data_generate as fdg
from methylbert.utils import set_seed
from torch.utils.data import DataLoader
from methylbert.data.vocab import MethylVocab
from methylbert.data.dataset import MethylBertFinetuneDataset
from methylbert.trainer import MethylBertFinetuneTrainer
import os
import pandas as pd
from methylbert.deconvolute import deconvolute
from transformers import AutoModel

''' This part preprosses the tumor and normal data later used for fine tuning the model.
    The output files are seq_train.csv, seq_test.csv and dmrs.csv.
    They can later be loaded into a dataloader.

'''

print("Fine Tuning start...")

f_bam_file_list = "../fine_tune_data.txt"
f_bam = "../data/bam_for_classification/sorted_bulk_data.bam"
f_dmr = "../data/reference/dmr_filtered.csv"
f_ref = "../data/reference/hg38.fa"
out_dir = "test/"


fdg.finetune_data_generate(
    sc_dataset = f_bam_file_list,
    f_dmr = f_dmr,
    f_ref = f_ref,
    output_dir=out_dir,
    split_ratio = 0.8, # Split ratio to make training and validation data
    n_mers=3, # 3-mer DNA sequences
    n_cores=20
)

# This part is for fine tuning the model.

set_seed(42)
seq_len=100
n_mers=3
batch_size=5
num_workers=0
output_path = "tmp/fine_tune/"

print("[1/13] set fine-tuning parameters")

# Creat a look-up table
tokenizer = MethylVocab(n_mers)

print("[2/13] Created Tokenizer")

# Load the data files int a data set object
train_dataset = MethylBertFinetuneDataset("tmp/train_seq.csv",
                                          tokenizer,
                                          seq_len=seq_len)
test_dataset = MethylBertFinetuneDataset("tmp/test_seq.csv",
                                         tokenizer,seq_len=seq_len)

print("[3/13] Data loaded to Dataset object.")

# Load the data into a data loader
train_data_loader = DataLoader(train_dataset, batch_size=batch_size, num_workers=num_workers, pin_memory=False, shuffle=True)

test_data_loader = DataLoader(test_dataset, batch_size=batch_size, num_workers=num_workers, pin_memory=True, shuffle=False)

print("[4/13] Data loaded into data loader.")

# creating trainer object for fine tuning model
trainer = MethylBertFinetuneTrainer(len(tokenizer),
                      save_path=output_path,
                      train_dataloader=train_data_loader,
                      test_dataloader=test_data_loader,
                      lr=1e-4, with_cuda=False,
                      log_freq=1,
                      #eval_freq=10, #activate this only when you want to evaluate the model with test_data_loader
                      warmup_step=3,
                      use_multiprocessing=False
                      )

print("[5/13] Created trainer object")

# load pretrained model into trainer object
trainer.load("hanyangii/methylbert_hg19_2l") # alternative numbers of encoder blocks: 2,4,6,8,12

print("[6/13] Loaded model into trainer.")

''' Training produces the following files in the tmp/fine_tune/ directory.

    config.json and pytorch_model.bin: model configuration and the trained MethylBERT model
    dmr_encoder.pickle : The trained DMR encoder in the MethylBERT model
    read_classification_model.pickle : The trained fully connected neural network for read classification
    train.csv and eval.csv : tracked training and evaluation loss and accuracy values during the training

'''
trainer.train(steps=1)
print("[7/13] Completed Training")


''' Second part of the MethylBERT procedure. Fine Tuning needs to be finished. '''

''' This part preprocesses the input bulk data used for testing later.
    Here, the output file data.csv is generated.
    This should also be loaded into a dataloader later.

'''

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


print("[8/13] Setup for deconvolution complete")


print("[9/13] Created tokenizer")

# loading the bulk data
dataset = MethylBertFinetuneDataset("tmp/data.csv",
                                    tokenizer,
                                    seq_len=seq_len)
data_loader = DataLoader(dataset, batch_size=batch_size, num_workers=num_workers)


print("[10/13] Created dataset")

restore_dir = "tmp/fine_tune/"
trainer = MethylBertFinetuneTrainer(len(tokenizer),
                                    train_dataloader=train_data_loader,
                                    test_dataloader=data_loader, with_cude=False
                                    )
print("[11/13] Trainer created")

#trainer.load("model") # Load the fine-tuned MethylBERT model
trainer.load(restore_dir)

print("[12/13] model loaded")

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

print("[13/13] Deconvolution complete.")
