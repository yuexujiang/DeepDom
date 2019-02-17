import sys
import os
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
from keras.layers import Activation
import random
from keras.optimizers import Adam
from keras.layers import Concatenate
from keras.models import load_model
from keras import layers
import keras.backend as K
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

def convertlabels_to_categorical(seq):
    
    letterDict = {} 
    letterDict["0"] = [1,0,0]
    letterDict["1"] = [0,1,0]
    letterDict["-"] = [0,0,1]
    AACategoryLen = 3 # 3 for '-'
    
    probMatr = np.zeros((len(seq),AACategoryLen))
    
    AANo  = 0
    for AA in seq:
        if not AA in letterDict:
           probMatr[AANo] = np.full(AACategoryLen, 0)
        elif AA=="-":
           probMatr[AANo]= letterDict[AA]
        else:
           probMatr[AANo]= [1-float(AA),float(AA),0]
        
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

def process_inputlabels(file):
    prot_id = ''
    prot_seq = ''
    labels = []
    ids = []
    for line in open(file, 'r'):
        if line[0] == '>':
            if prot_id != '':
                labels.append(prot_seq)
                ids.append(prot_id)
            
            prot_id = line.strip()
            prot_seq = ''
        elif line.strip() != '':
                prot_seq = prot_seq + line.strip()
    
    if prot_id != '':
        labels.append(prot_seq)
        ids.append(prot_id)
    return (ids,labels)

def main():
    parser=argparse.ArgumentParser(description='MusiteDeep prediction tool for general, kinase-specific phosphorylation prediction or custom PTM prediction by using custom models.')
    parser.add_argument('-seqfile',  dest='seqfile', type=str, help='Processed sequence data to be trained (output from dataprocess.pl).', required=True)
    parser.add_argument('-labelfile',  dest='labelfile', type=str, help='Processed label data to be trained (output from dataprocess.pl).', required=True)
    parser.add_argument('-model-prefix',  dest='modelprefix', type=str, help='File name of the generated model.', required=True)
    parser.add_argument('-lossfile',  dest='lossfile', type=str, help='Specify a file to store "loss" history during the training process', required=False,default="history_loss.txt")
    args = parser.parse_args()
    seqfile=args.seqfile; #file name of the processed sequence data (output from dataprocess.pl)
    labelfile=args.labelfile; #file name of the processed label data (output from dataprocess.pl)
    modelprefix=args.modelprefix
    lossfile = args.lossfile
    (ids,seqs)=process_inputseqs(seqfile)
    (ids,labels)=process_inputlabels(labelfile)
    rawdata=zip(seqs,labels)
    random.shuffle(rawdata)
    inputX=[convertSampleToPhysicsVector_pca(i[0]) for i in rawdata]
    inputY=[convertlabels_to_categorical(i[1]) for i in rawdata]
    
    
    data=zip(inputX,inputY)
    random.seed(4)
    random.shuffle(data)
    random.shuffle(data)
    
    train_num=int(len(inputX)*0.9)
    train=data[0:train_num]
    val=data[train_num:]
    
    
    
    trainX=[i[0] for i in train]
    trainY=[i[1] for i in train]
    xx=np.dstack(trainX)
    xx=np.rollaxis(xx,-1)
    
    yy=np.dstack(trainY)
    yy=np.rollaxis(yy,-1)
    
    valX=[i[0] for i in val]
    valY=[i[1] for i in val]
    valX=np.dstack(valX)
    valX=np.rollaxis(valX,-1)
    valY=np.dstack(valY)
    valY=np.rollaxis(valY,-1)
    
    ####################
    input =Input(shape=(None,inputX[0].shape[1]))
    
    x1 = layers.Bidirectional(CuDNNLSTM(100,return_sequences=True),merge_mode='sum')(input)
    
    x2 = layers.Bidirectional(CuDNNLSTM(100,return_sequences=True),merge_mode='sum')(x1)
    
    x3 = layers.Bidirectional(CuDNNLSTM(100,return_sequences=True),merge_mode='sum')(x2)
    x4 = layers.Bidirectional(CuDNNLSTM(3,return_sequences=True),merge_mode='sum')(x3)
    
    output=Activation('softmax')(x4)
    
    model=Model(input,output)
    early_stopping = EarlyStopping(monitor='val_loss', patience=10)
    model.compile(loss='binary_crossentropy',optimizer='adam',metrics=['accuracy'])
    
    vv=(valX,valY)
    
    fitHistory_batch = model.fit(x=xx,y=yy, batch_size=1000, validation_data=vv,epochs=100,callbacks=[early_stopping])
    
    
    model.save(modelprefix)                                           #specify the location where your trained model(.h5 file) would be saved
    
    trainloss=fitHistory_batch.history['loss']
    valloss=fitHistory_batch.history['val_loss']
    outputfile = open(lossfile, 'w')                 #specify a file to store "loss" history during the training process
    
    for i in range(len(trainloss)):
      print >> outputfile, "%f\t%f"%(trainloss[i],valloss[i])
    
    outputfile.close()
    
    

if __name__ == "__main__":
    main()         
