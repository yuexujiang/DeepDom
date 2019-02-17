# DeepDom
DeepDom is an ab-initio method for protein domain boundary prediction

Usage:

1. Train a model by your own data

To train a model, users should have .fasta file for protein sequences and domain boundary annotations on corresponding protein sequences. Then use the "dataprocess.pl" code to transformat the data. After that, the processed data can be used as input for "train.py" to train a model. Note: the requirement of packages that imported in the code need to be met.

2. Use the existing trained model that mentioned in our paper

By fine tuning the hyper-parameters, We have already trained a model and saved it as file "foldmodel_bilstmwrapper_4sum200_80_40nr_sliwin.h5".

3. Predict

To predict domain boundary for protein sequences, firstly, users need to transformat the .fasta sequence using "dataprocess.pl" (run "perl dataprocess.pl -h" to see the helps) and using "predict.py" to predict for protein sequences (run "python predict.py -h" to see the helps). Either users' own model or the model we provided can be used for prediction.
Examples:

perl dataprocess.pl -input_seq sample_data_seq.txt -input_label sample_data_label.txt -output_seq processed_seq.txt -output_label processed_label.txt

python predict.py -input processed_seq.txt -output predict_result.txt


If you find the codes or our method is useful, please cite our paper "DeepDom: Predicting protein domain boundary from sequence alone using stacked bidirectional LSTM".
