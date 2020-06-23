###########################################
#☄SARS-CoV-2 sequence data processing
# Author: Dong Liang
# Email: ldifer@gmail.com
# May 21, 2020
###########################################


# import pysam
# import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
import tensorflow as tf 
from tensorflow.keras.preprocessing.sequence import pad_sequences 
import numpy as np
import pandas as pd
from collections import defaultdict
import re
import yaml


def load_virus_info(type = 'nucleotide', filter_on = None, data_preprocessing = True, col = None):

    path = './' + type + '.csv'
    
    # Load the dataset   
    virus_info = pd.read_csv(path, index_col = False)
    
    if data_preprocessing:
        # Set index    
        virus_info.reset_index(drop = True, inplace = True)

        # Data proprocessing
        virus_info['Collection_Date'] = pd.to_datetime(virus_info['Collection_Date'])

        # Locations
        virus_info['Country_Region'] = virus_info['Geo_Location'].str.split(':').str.get(0).str.upper().str.strip()
        virus_info['City'] = virus_info['Geo_Location'].str.split(':').str.get(1).str.upper().str.strip()
        virus_info['City'].replace(to_replace = {
            'MICHIGAN': 'MI', 
            'ILLINOIS': 'IL',
            'WISCONSIN': 'WI', 
            'OREGON, WISCONSIN': 'WI', 
            'DANE COUNTY, WISCONSIN': 'WI', 
            'MOUNT HOREB, WISCONSIN': 'WI', 
            'MADISON, WISCONSIN': 'WI', 
            'VERONA, WISCONSIN': 'WI', 
            'MADISON, WISCONSIN': 'WI', 
            
            'MARINGOUIN, LA': 'LA', 
            'SAINT ROSE, LA': 'LA', 
            'NEW ORLEANS, LA': 'LA', 
            'LOCKPORT, LA': 'LA', 
            'LULING, LA': 'LA', 
            'RACELAND, LA': 'LA', 
            
            'KENNER, LA': 'LA', 
            'GHEENS, LA': 'LA', 
            'THIBODAUX,LA': 'LA', 
            'LACOMBE, LA': 'LA', 
            'HOUMA, LA': 'LA', 
            'SLIDELL, LA': 'LA', 
            'EAST FELICIANA PARISH, LOUISIANA': 'LA',
            
            'NEW YORK': 'NY', 
            'ARIZONA': 'AZ', 
            'NORTH CAROLINA': 'NC', 
            'SAN FRANCISCO, CA': 'CA', 
            'CA, SAN DIEGO COUNTY': 'CA', 
            'SNOHOMISH COUNTY, WA': 'WA'
            
        }, inplace = True)

    # Columns
    if isinstance(col, list):
        virus_info = virus_info.loc[:, col]

    return virus_info


def load_virus_seq(path):
    # Initialization
    records = SeqIO.parse(path, 'fasta')    
    seqs = {}
    for record in records:
        id = record.id.split('.')[0]
        if id.find('join') == -1:
            seqs[id] = {}

    # Extract nucleotide sequences from open reading frames 
    records = SeqIO.parse(path, 'fasta')    
    for record in records:
        id = record.id.split('.')[0]
        seq = str(record.seq)
        ORF = record.description.split('|')[1].split(' [')[0]
        if id.find('join') == -1:
            seqs[id].update({ORF: seq})
        
    return seqs



class COVID19(object):
    def __init__(self, id, seq, is_coding = True):
        self.code = 'ATCG'
        self.id = id
        self.seq = seq  # {}
        self.is_coding = is_coding
        self.seq_spike = self.get_S_seq()
        if self.is_coding:
            pass
        

    def one_hot_encoding(self, seq):
        # Version 1
        mapping = dict(zip(self.code, range(4)))
        seq_ohe = [mapping.get(s, 4) for s in seq]
        return np.eye(5)[seq_ohe]

        # Version 2
        # seq = list(seq)
        # label_encoder = LabelEncoder()
        # index_seq = label_encoder.fit_transform(list(seq))
        # ohe = OneHotEncoder(handle_unknown='ignore')
        # ohe_seq = ohe.fit_transform(index_seq.reshape(-1, 1))
        # return ohe_seq.toarray()

    def normalization(self):
        pass

    def get_S_seq(self):

        # Extract S protein sequences
        if not self.is_coding: 
            spike_start = 'atgtttgttt'.upper()
            spike_end = 'cattacacataa'.upper()

            left_ind = self.seq.find(spike_start)
            right_ind = self.seq.find(spike_end) + len(spike_end)

            spike_seq = self.seq[left_ind:right_ind]
        else:
            if 'surface glycoprotein' in self.seq.keys():
                spike_seq = self.seq['surface glycoprotein']
            else: 
                # print('not present')
                spike_seq = ''
        
        return spike_seq
        # assert len(spike_seq) == 3822

    def get_SNP_genotype(self):
        pass

    def get_dN_dS(self):
        pass
    


    def __call__(self, ORF = 'spike', ohe = True, return_residual = False):

        if ORF == 'complete':
            seq = list(self.seq.values())[0]
            seq_ohe = self.one_hot_encoding(seq)
            return seq_ohe

        elif ORF == 'spike':
            # One hot encoding
            seq_spike_ohe = self.one_hot_encoding(self.seq_spike)
            self.seq_spike_ohe = seq_spike_ohe
            
            if return_residual:
                # Reference sequences for Spike 
                a_yaml_file = open("./covid19_toolkit/reference_seq.yaml")
                parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)
                seq_spike_ref = parsed_yaml_file['spike_protein'].strip()
                # print(len(seq_spike_ref))
                seq_spike_ref_ohe = self.one_hot_encoding(seq_spike_ref)

                # Residual sequence
                try:
                    seq_spike_ohe_resid = self.seq_spike_ohe - seq_spike_ref_ohe
                except:
                    seq_spike_ohe_resid = None
                
                self.seq_spike_ohe_resid = seq_spike_ohe_resid
                
                return seq_spike_ohe_resid
            else:
                return seq_spike_ohe 



def preprocessing(dataset):
    covid19_seqs = [COVID19(id, seq) for id, seq in dataset.items()]
    covid19_complete = {seq.id: seq(ORF = 'complete', return_residual = False) for seq in covid19_seqs }
    covid19_dataset = zero_padding(covid19_complete)
    return covid19_dataset

def zero_padding(dataset):
    # Zero padding 
    maxlen = max([d.shape[0] for k, d in dataset.items()])
    dataset_padded = {k: pad_sequences(d.T, maxlen = maxlen, padding = 'post', dtype = 'int').T for k, d in dataset.items()}
    return dataset_padded
    

class Feature_engineering():
    def extract_SNP(self):
        pass

    def get_S_protein(self):
        pass

    def get_dN_dS(self):
        pass
