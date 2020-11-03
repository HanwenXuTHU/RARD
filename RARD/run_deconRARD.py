from RARD.DeconRARD import DeconRARD
import pandas as pd
import collections
import numpy as np


def deconvolution(ref_path, mix_path, marker_path, save_path, unknown=False, scale=1):
    """
    deconvolution by DeconCCN
    :param ref_path: path to ref.csv, format should be same as the demo
    :param mix_path: path to mix.csv, format should be same as the demo
    :param marker_path: path to markers.csv, format should be same as the demo
    :param save_path: path to save the deconvolution results
    :param unknown: unknown content
    :return:
    """
    ref = pd.read_csv(ref_path, index_col=0)
    mix = pd.read_csv(mix_path, index_col=0)
    markers = pd.read_csv(marker_path, index_col=0)
    markers = markers.index.values
    cell_type = ref.columns.values
    samples = mix.columns.values
    prop = collections.OrderedDict()
    prop['cell types'] = cell_type
    reference = []
    mixture = []
    for i in range(len(markers)):
        reference.append(ref.loc[markers[i]])
        mixture.append(mix.loc[markers[i]])
    reference = np.asarray(reference)
    mixture = np.asarray(mixture)
    prop_predict = DeconRARD(scale * reference, scale * mixture, unknown=unknown)
    for i in range(len(samples)):
        prop[samples[i]] = []
        for j in range(len(cell_type)):
            prop[samples[i]].append(prop_predict[j, i])
    prop = pd.DataFrame(prop)
    prop.to_csv(save_path, index=False)

