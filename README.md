# RARED: Robust and accurate rare cell type de-convolution for gene expression and DNA methylation data

![Image text](https://github.com/HanwenXuTHU/RARD/tree/master/img/decon.png)

## Section 1: System Requirements
&emsp; RARD was implemented using python. It can be installed in Windows, Linux as well as MAC OS. RARD requires python version >= 3 and all the dependent packages will be installed using pip.

## Section 2: Installation Instruction

&emsp; RARD can be installed using pip by the following command:

``` shell
pip install RARD
```


## Section 3: How to Use RARD

### Section 3.1: Input Data Preparation

&emsp; We provide demon [dataset](https://github.com/HanwenXuTHU/RARD/tree/master/demo_data) generated from single cell mouse RNA-seq.

&emsp; Note: RARD is not sensitive to gene identifier, but the user should restrict gene identifier to a specific one.


### Section 3.2: Deconvolution

We provide a demo here to illustrate how to run RARD

``` Python
from RARD.run_deconRARD import deconvolution

deconvolution('ref.csv', 'mix.csv', 'marker.csv', 'prop_predict.csv', scale=0.01)
```

The results will be saved in prop_predict.csv.


