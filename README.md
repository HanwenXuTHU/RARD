# Aida: accurate inference of cell type fraction across various data sources

## Section 1: System Requirements
&emsp; Aida was implemented using python. It can be installed in Windows, Linux as well as MAC OS. DeconCCN requires python version >= 3 and all the dependent packages will be installed using pip.

## Section 2: Installation Instruction

&emsp; DeconCCN can be installed using pip by the following command:

``` shell
pip install AidaGo
```


## Section 3: How to Use Aida

### Section 3.1: Input Data Preparation

&emsp; We provide demon [dataset](https://github.com/HanwenXuTHU/Aida/tree/master/demo_data) generated from single cell mouse RNA-seq.

&emsp; Note: DeconCCN is not sensitive to gene identifier, but the user should restrict gene identifier to a specific one.


### Section 3.2: Deconvolution

We provide a demo here to illustrate how to run Aida

``` Python
from Aida.run_deconccn import deconvolution

deconvolution('ref.csv', 'mix.csv', 'marker.csv', 'prop_predict.csv', scale=0.01)
```

The results will be saved in prop_predict.csv.


