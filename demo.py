# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:18 2019

@author: Hanwen Xu

E-mail: hw-xu16@mails.tsinghua.edu.cn
"""
import numpy as np
from DeconCCNMethod import DeconCCN
import matplotlib.pyplot as plt

reference = np.load("reference.npy")
demoMixture = np.load("demoMixture.npy")
proportionDeconvolution = DeconCCN(reference, demoMixture)
cellType = ['Neutrophils', 'NK', 'B cells', 'CD4+', 'Monocytes', 'CD8+']
for i in range(6):
    rects = plt.bar(1 + i * 0.2, proportionDeconvolution[i, 0], width = 0.05)
plt.xticks([1.0, 1.2, 1.4, 1.6, 1.8, 2.0], cellType)
plt.show()
debug = 0