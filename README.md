# DeconCCN:a cell subtype deconvolution algorithm based on component-wise condition number

## 1.System requirements

We developed DeconCCN on Windows 10 platform. DeconCCN requires python3 and all the packages will be installed using pip.

## 2.Installation guide

DeconCCN should be installed using pip by following the command below:

`pip install DeconCCN`

This may take a few minutes, which depends on the required packages.

## 3.Demo

### 3.1 data

Data used in the demo is the dataset we test as single cell mouse dataset.

Data has been uploaded to the [demo_data](https://github.com/HanwenXuTHU/DeconCCN/tree/master/demo_data) file.

### 3.2 run the demo

We provide a demo here to illustrate how to run 

`from DeconCCN.run_deconccn import deconvolution`

`deconvolution('ref.csv', 'mix.csv', 'marker.csv', 'prop_predict.csv')`

The results will be saved in prop_predict.csv.
