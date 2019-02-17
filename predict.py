import sys
import re
import pandas as pd
import keras.utils.np_utils as kutils
from keras.layers import GlobalMaxPooling1D
from keras.layers import Input, Dense, Conv1D, MaxPooling1D, UpSampling1D,GlobalMaxPooling1D,Dropout,BatchNormalization
from keras.layers import *
from keras.layers import CuDNNLSTM
import keras
from keras.models import Model, Sequential
from  keras.callbacks import EarlyStopping, ModelCheckpoint
import numpy as np
#from keras.activations import softmax
from keras.layers import Activation
import random
from keras.optimizers import Adam
from keras.layers import Concatenate
from keras.models import load_model
from keras import layers
import keras.backend as K
from numpy import median
import argparse
def convertSampleToPhysicsVector_pca(seq):
    """
    Convertd the raw data to physico-chemical property
    
    PARAMETER
    seq: "MLHRPVVKEGEWVQAGDLLSDCASSIGGEFSIGQ" one fasta seq
        X denoted the unknow amino acid.
    
    
    probMatr: Probability Matrix for Samples. Shape (nb_samples, 1, nb_length_of_sequence, nb_AA)
    """
    
    letterDict = {} 
    letterDict["A"] = [0.008,0.134,-0.475,-0.039,0.181,0]
    letterDict["R"] = [0.171,-0.361,0.107,-0.258,-0.364,0]
    letterDict["N"] = [0.255,0.038,0.117,0.118,-0.055,0]
    letterDict["D"] = [0.303,-0.057,-0.014,0.225,0.156,0]
    letterDict["C"] = [-0.132,0.174,0.070,0.565,-0.374,0]
    letterDict["Q"] = [0.149,-0.184,-0.030,0.035,-0.112,0]
    letterDict["E"] = [0.221,-0.280,-0.315,0.157,0.303,0]
    letterDict["G"] = [0.218,0.562,-0.024,0.018,0.106,0]
    letterDict["H"] = [0.023,-0.177,0.041,0.280,-0.021,0]
    letterDict["I"] = [-0.353,0.071,-0.088,-0.195,-0.107,0]
    letterDict["L"] = [-0.267,0.018,-0.265,-0.274,0.206,0]
    letterDict["K"] = [0.243,-0.339,-0.044,-0.325,-0.027,0]
    letterDict["M"] = [-0.239,-0.141,-0.155,0.321,0.077,0]
    letterDict["F"] = [-0.329,-0.023,0.072,-0.002,0.208,0]
    letterDict["P"] = [0.173,0.286,0.407,-0.215,0.384,0]
    letterDict["S"] = [0.199,0.238,-0.015,-0.068,-0.196,0]
    letterDict["T"] = [0.068,0.147,-0.015,-0.132,-0.274,0]
    letterDict["W"] = [-0.296,-0.186,0.389,0.083,0.297,0]
    letterDict["Y"] = [-0.141,-0.057,0.425,-0.096,-0.091,0]
    letterDict["V"] = [-0.274,0.136,-0.187,-0.196,-0.299,0]
    letterDict["X"] = [0,-0.00005,0.00005,0.0001,-0.0001,0]
    letterDict["-"] = [0,0,0,0,0,1]
    AACategoryLen = 6 #6 for '-'
    probMatr = np.zeros((len(seq),AACategoryLen))
    AANo  = 0
    for AA in seq:
        if not AA in letterDict:
           probMatr[AANo] = np.full(AACategoryLen, 0)
        else:
           probMatr[AANo]= letterDict[AA]
        
        AANo += 1
    
    return probMatr

def process_inputseqs(file):
    prot_id = ''
    prot_seq = ''
    seqs = []
    ids = []
    for line in open(file, 'r'):
        if line[0] == '>':
            if prot_id != '':
                seqs.append(prot_seq)
                ids.append(prot_id)
            
            prot_id = line.strip()
            prot_seq = ''
        elif line.strip() != '':
                prot_seq = prot_seq + line.strip()
    
    if prot_id != '':
        seqs.append(prot_seq)
        ids.append(prot_id)
    return (ids,seqs)


def initial_list(length):
    target=[None]*length
    for ind in range(len(target)):
        target[ind]=[]
    return target
    


def update_list(stride,num,comp,pred_score,pred):
    for j in range(comp):
        pred[j+num*stride].append(pred_score[j])  
    return pred  



def main():
    parser=argparse.ArgumentParser(description='MusiteDeep prediction tool for general, kinase-specific phosphorylation prediction or custom PTM prediction by using custom models.')
    parser.add_argument('-input',  dest='inputfile', type=str, help='Processed sequence data to be predicted (output from dataprocess.pl).', required=True)
    parser.add_argument('-output',  dest='outputfile', type=str, help='file name of the prediction results.', required=True)
    parser.add_argument('-model-prefix',  dest='modelprefix', type=str, help='File name of the custom model used for prediciton. If donnot have one, please run train.py to generate the model', required=False,default=None)
    parser.add_argument('-windowsize', dest='window',type = int, help='specify the window size, must be the same values as in dataprocess.pl',required=False,default=300)
    parser.add_argument('-stridesize', dest='stride',type = int, help='specify the stride size, must be the same values as in dataprocess.pl',required=False,default=80)
    args = parser.parse_args()
    inputfile=args.inputfile;
    outputfile=args.outputfile;
    modelprefix=args.modelprefix
    win = args.window
    stride = args.stride
    ##########################To predict#############################
    if modelprefix is None:
         modelprefix = "foldmodel_bilstmwrapper_4sum200_80_40nr_sliwin.h5"
    
    model=load_model(modelprefix)
    (testids,testseqs)=process_inputseqs(inputfile)
    testX=[convertSampleToPhysicsVector_pca(i) for i in testseqs]
    
    tt=np.dstack(testX)
    tt=np.rollaxis(tt,-1)
    posscore=model.predict(tt)
    fo=open(outputfile, 'w')
    
    regex = r"(\d+)left(\d+)"
    pred=""
    for i in range(len(testids)):
         com=testids[i].split("_")
         seql=int(com[3])
         num_ind=com[4]  
         pred_score=posscore[i][:,1]  
         if re.search(regex,num_ind):
             match = re.search(regex, num_ind)
             num=int(match.group(1))
             comp=int(match.group(2))   
             if num==0:
               pred=np.zeros(seql)
               name="_".join(p for p in com[0:4])
               fo.write(name+"\n")
             for j in range(comp):
                 if pred_score[j]>pred[j+num*stride]:
                      pred[j+num*stride]=pred_score[j]   
             score=' '.join([str(x) for x in pred])   
             fo.write(score+"\n")
         elif num_ind=="0":
             pred=np.zeros(seql)
             num=int(num_ind)
             for j in range(100):
                 if pred_score[j]>pred[j+num*stride]:
                      pred[j+num*stride]=pred_score[j]
             name="_".join(p for p in com[0:4])
             fo.write(name+"\n")    
             if (num*stride+win)==seql:
                 score=' '.join([str(x) for x in pred])
                 fo.write(score+"\n")
         else:    
             num=int(num_ind)
             for j in range(100):
                 if pred_score[j]>pred[j+num*stride]:
                      pred[j+num*stride]=pred_score[j]
             if (num*stride+win)==seql:
                 score=' '.join([str(x) for x in pred])
                 fo.write(score+"\n")
    
    
             
    fo.close() 


if __name__ == "__main__":
    main()         
