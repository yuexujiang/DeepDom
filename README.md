# DeepDom
DeepDom is an ab-initio method for protein domain boundary prediction

Usage:

1. Train a model by your own data

To train a model, users should have .fasta file for protein sequences and domain boundary annotations on corresponding protein sequences. Then use the "dataprocess.pl" code to transformat the data. After that, the processed data can be used as input for "train.py" to train a model. Note: the requirement of packages that imported in the code need to be met.

2. Use the existing trained model that mentioned in our paper

By fine tuning the hyper-parameters, We have already trained a model and saved it as file "foldmodel_bilstmwrapper_4sum200_80_40nr_sliwin.h5".

3. Predict

To predict domain boundary for a protein sequence, firstly, users need to transformat the .fasta sequence using "dataprocess.pl" as usage 1 and specify the location of the processed sequence data in "predict.py". Users also need to load a trained model, either users' own model or the model we provided.

We comment the parts in the codes where users' inputs are needed.

If you find the codes or our method is useful, please cite our paper "DeepDom: Predicting protein domain boundary from sequence alone using stacked bidirectional LSTM".
