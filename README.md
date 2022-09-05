# DeepDom
DeepDom is an ab-initio method for protein domain boundary prediction
## Installation

  - Installation has been tested in Linux and Mac OS X with Python 3.5. 
  - Since the package is written in python 3.5, [python 3.5](https://www.python.org/downloads/ ) with the pip tool must be installed first. 
DeepDom uses the following dependencies:
numpy,  scipy, pandas, h5py, keras version==2.3.0, tensorflow==1.13.1
You can install these packages first, by the following commands:

```sh
pip install pandas
pip install numpy
pip install scipy
pip install h5py
pip install -v keras==2.1.2
pip install tensorflow (or GPU supported tensorflow, refer to https://www.tensorflow.org/install/ for instructions)
```
 - This is the Tensorflow version, you must change the backend to TensorFlow.
If you have run Keras at least once, you will find the Keras configuration file at:
$HOME/.keras/keras.json
If it isn’t there, you can create it. 
Change the default configuration file into:
```sh
{	
    "image_dim_ordering": "th",
    "epsilon": 1e-07,
    "floatx": "float32",
    "backend": "tensorflow"
}
```
## Running on GPU or CPU

>If you want to use GPU, you also need to install [CUDA]( https://developer.nvidia.com/cuda-toolkit) and [cuDNN](https://developer.nvidia.com/cudnn); refer to their websites for instructions. 
CPU is only suitable for prediction not training. 
## Training data

Since the data size is larger than the Github limitation, the training data for DeepDom is stored at Google Drive that can be accessed through the link https://drive.google.com/drive/folders/1HsHEbk1gjtEdivr9vvOpXPK7MsoBZp3h?usp=sharing
## Usage:
## 1. Train a model by your own data

To train a model, users should have .fasta file for protein sequences and domain boundary annotations on corresponding protein sequences.
Please refer to the sample_data_seq.txt and sample_data_label.txt as examples. Then use the "dataprocess.pl" to transformat the data. (run "perl dataprocess.pl -h" to see the helps) After that, the processed data can be used as input for "train.py" to train a model (run "python train.py -h" to see the helps). Note: the requirement of packages that imported in the code need to be met. 
#### Examples (using our provided example data): 
 
```sh
perl dataprocess.pl -input_seq sample_data_seq.txt -input_label sample_data_label.txt -output_seq processed_seq.txt -output_label processed_label.txt

python train.py -seqfil processed_seq.txt -labelfile processed_label.txt -model-prefix custom_model.h5
```

custom_model.h5 is the model generated, users can use this file to predict and can also use our pre-trained model that mentioned in our paper. The pre-trained model was saved in file "foldmodel_bilstmwrapper_4sum200_80_40nr_sliwin.h5".

## 2. Predict

To predict domain boundary for protein sequences, firstly, users need to transformat the .fasta sequence using "dataprocess.pl" (run "perl dataprocess.pl -h" to see the helps) and using "predict.py" to predict for protein sequences (run "python predict.py -h" to see the helps). Either users' own model or the model we provided can be used for prediction.
#### Examples (using our provided example data):
For GPU users, predict.py will use our pre-trained GPU model "foldmodel_bilstmwrapper_4sum200_80_40nr_sliwin.h5":
```sh
perl dataprocess.pl -input_seq sample_data_seq.txt -output_seq processed_seq.txt

python predict.py -input processed_seq.txt -output predict_result.txt
```
For CPU users, be sure to add -model-prefix cpu_model.h5 when call predict.py:
```sh
perl dataprocess.pl -input_seq sample_data_seq.txt -output_seq processed_seq.txt

python predict.py -input processed_seq.txt -output predict_result.txt -model-prefix cpu_model.h5
```
Or users can use the custom model trained (as shown in 1) by their own data to predict.
 ```sh
python predict.py -input processed_seq.txt -output predict_result.txt -model-prefix custom_model.h5
```

========================================================================================================================================
#### If you find the codes or our method is useful, please cite our paper "DeepDom: Predicting protein domain boundary from sequence alone using stacked bidirectional LSTM".
