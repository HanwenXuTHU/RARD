# -*- coding: utf-8 -*-
"""
Created on Wen Apr 8 10:18 2020

@author: Hanwen Xu

E-mail: hw-xu16@mails.tsinghua.edu.cn

New update: process data with unknown content
"""
import numpy as np
from tqdm import tqdm
from scipy import optimize
from sklearn.svm import NuSVR
import statsmodels.api as sm


def DeconCCN(reference0, markerData0, unknown=False):
    """
    DeconCCN: a cell subtype deconvolution algo-rithm based on component-wise condition number
    :param reference0: an N×K matrix, which denotes the external reference panel, should be ndarray.
    :param markerData0: N × M matrix, which denotes the mixture data, should be ndarray
    :param unknown: whether or not contain the unknown content
    :return: results of deconvolution
    """
    reference,markerData = filt_zeros(reference0, markerData0)
    markerNum = np.size(reference, 0)
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(markerData, 1)
    bestReference = []
    bestMarker = []
    conditionNumberHistory = []
    bestNumberHistory = []
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    for i in tqdm(range(numOfSamples)):
        selectArea = np.arange(markerNum)
        selectMixture = markerData[selectArea, i]
        selectReference = reference[selectArea, :]
        minimumConditionNumber = 10 ** (100)
        endNumber = np.size(selectReference, 0)
        for selectedNumber in range(int(endNumber / 10)):
            minDistances = 10 ** (50)
            for j in range(tissueNum):
                for k in range(j + 1, tissueNum):
                    distances = selectReference[:, j] - selectReference[:, k]
                    distances = np.sqrt(np.sum(np.multiply(distances, distances)))
                    if distances < minDistances:
                        minDistances = distances
                        closetJ = j
                        closetK = k
            sumData = selectReference[:, closetJ] + selectReference[:, closetK]
            area = sumData == 0
            sumData[area] = 10**(-100)
            collinearity = np.abs(selectReference[:, closetJ] - selectReference[:, closetK]) / (sumData)
            collinearityIndex = np.argsort(collinearity)
            area = np.ones(np.size(selectReference, 0))
            area = area.astype(np.bool)
            area[collinearityIndex[0]] = False
            selectArea = np.arange(np.size(selectReference, 0))
            selectArea = selectArea[area]
            selectMixture = selectMixture[selectArea]
            selectReference = selectReference[selectArea, :]
            ConditionNumber = mixedConditionNumber(selectReference, selectMixture)
            conditionNumberHistory.append(ConditionNumber)
            if ConditionNumber < minimumConditionNumber:
                minimumConditionNumber = ConditionNumber
                bestReference = selectReference
                bestMarker = np.zeros([np.size(selectReference, 0), 1])
                bestMarker[:, 0] = selectMixture
                bestNumber = selectedNumber
        t = RobustSVR(bestReference, bestMarker, unknown=unknown)
        bestNumberHistory.append(bestNumber)
        proportionDeconvolution[:, i] = t[:, 0]
    return proportionDeconvolution


def mixedConditionNumber(referenceSelect, markerSelect):
    """
    caculate the largest component-wise condition number
    :param referenceSelect: reference data in selected markers
    :param markerSelect: mixture data in selected markers
    :return: the largest component-wise condition number.
    """
    pinvReferenceSelect = np.linalg.pinv(referenceSelect)
    bNorm = np.linalg.norm(markerSelect)
    maxConditionNumber = 0
    tissueNumber = np.size(referenceSelect, 1)
    for i in range(tissueNumber):
        tq = pinvReferenceSelect[i, :]
        conditionNumber = (bNorm / np.abs(np.dot(tq, markerSelect))) * np.linalg.norm(tq)
        if conditionNumber > maxConditionNumber:
            maxConditionNumber = conditionNumber
    return maxConditionNumber


def RobustSVR(reference, markerData, unknown=False):
    """
    outlier detection and support vector regression
    :param reference: reference after marker selection
    :param markerData: mixture after marker selection
    :param unknown: whether or not contain unknown content
    :return: results of deconvolution
    """
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(markerData, 1)
    markerNumber = np.size(reference, 0)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    for i in range(numOfSamples):
        iterNumber = 5
        iterReference = reference
        iterMarkerData = markerData[:, i]
        mixture = sm.RLM(iterMarkerData, iterReference).fit()
        test = mixture.params
        t = test / np.sum(test) if unknown == False else test
        c1 = np.zeros([np.size(iterReference, 0), 1])
        c1[:, 0] = iterMarkerData[:]
        t1 = np.zeros([tissueNum, 1])
        t1[:, 0] = t
        c2 = np.dot(iterReference, t1)
        res = c1 - c2
        s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
        res_std = np.abs(res / np.sqrt(s_2))
        res_Sort = np.sort(res_std[:, 0])
        T = res_Sort[int(0.75 * np.size(res_Sort))]
        memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
        memRef[:, :] = iterReference[:, :]
        memMix = np.zeros(np.size(iterMarkerData))
        memMix[:] = iterMarkerData[:]
        for j in range(iterNumber):
            mixture = sm.RLM(iterMarkerData, iterReference).fit()
            test = mixture.params
            t = test / np.sum(test) if unknown == False else test
            c1 = np.zeros([np.size(iterReference, 0), 1])
            c1[:, 0] = iterMarkerData[:]
            t1 = np.zeros([tissueNum, 1])
            t1[:, 0] = t
            c2 = np.dot(iterReference, t1)
            res = c1 - c2
            s_2 = np.sum(np.power(res, 2)) / (np.size(iterReference, 0) - np.size(iterReference, 1))
            res_std = res / np.sqrt(s_2)
            iterSelected = np.arange(np.size(iterReference, 0))
            area = np.abs(res_std[:, 0]) <= T
            iterSelected = iterSelected[area]
            iterReference = iterReference[iterSelected, :]
            iterMarkerData = iterMarkerData[iterSelected]
            if np.size(iterReference, 0) < int(tissueNum):
                iterReference = memRef
                iterMarkerData = memMix
                break
            if np.size(iterReference, 0) < int(0.5 * markerNumber):
                break
            memRef = np.zeros([np.size(iterReference, 0), np.size(iterReference, 1)])
            memRef[:, :] = iterReference[:, :]
            memMix = np.zeros(np.size(iterMarkerData))
            memMix[:] = iterMarkerData[:]
        weights = weightsDesigner(iterReference, iterMarkerData, unknown=unknown)
        weights = np.diag(weights)
        reference_w = np.dot(weights, iterReference)
        singleMarker_w = np.zeros([np.size(iterReference, 0), 1])
        singleMarker_w[:, 0] = np.dot(weights, iterMarkerData)
        t = nuSVR(reference_w, singleMarker_w, unknown=unknown)
        t = t[:, 0]
        if unknown == False:
            proportionDeconvolution[:, i] = t
        else:
            proportionDeconvolution[0 : -1, i] = t
            proportionDeconvolution[-1, i] = max(1 - np.sum(t), 0)
    return proportionDeconvolution


def weightsDesigner(ref, mix, unknown=False):
    """
    design a weights matrix to help us pay different attention on different markers.
    :param ref: reference
    :param mix: mixture
    :param unknown: unknown content
    :return: weights matrix
    """
    tissueNum = np.size(ref, 1)
    mixture = optimize.nnls(ref, mix)
    test = mixture[0]
    t = test / np.sum(test) if unknown == False else test
    c1 = np.zeros([np.size(ref, 0), 1])
    c1[:, 0] = mix[:]
    t1 = np.zeros([tissueNum, 1])
    t1[:, 0] = t
    small_cell = []
    for i in range(np.size(t1, 0)):
        if t1[i, 0] < 1:
            small_cell.append(i)
    if len(small_cell) > 0:
        scores = np.zeros([np.size(ref, 0), len(small_cell)])
        for j in range(len(small_cell)):
            for i in range(np.size(ref, 0)):
                temp = ref[i, :]
                temp1 = (temp - np.mean(temp)) / np.var(temp)
                scores[i, j] = np.abs(temp1[j]) / ((t1[j,0] + 0.001)**0.1)
        scores = np.max(scores, axis=1)
        scores = scores * 2 / (np.max(scores) - np.min(scores))
        scores = scores + (1 - np.min(scores)) * np.ones(np.size(scores))
        return scores
    else:
        return np.ones(np.size(ref, 0))


def nuSVR(reference, markerData, unknown=False):
    """
    nuSVR
    :param reference:
    :param markerData:
    :param unknown:
    :return:
    """
    nu = [0.25,0.50,0.75]
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(markerData, 1)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    nuHistory = []
    p0 = np.zeros([3, tissueNum, numOfSamples])
    for i in range(0, 3, 1):
        nuI = nu[i]
        clf = NuSVR(nu=nuI, kernel='linear')
        for j in range(numOfSamples):
            clf.fit(reference, markerData[:, j])
            t = clf.coef_
            t1 = np.zeros(tissueNum)
            t1[:] = t[0, :]
            area = t1 < 0
            t1[area] = 0
            t1 = t1 / np.sum(t1) if unknown == False else t1
            p0[i, :, j] = t1[:]
    for i in range(numOfSamples):
        minRMSE = 10**(50)
        truth = np.zeros([np.size(reference, 0), 1])
        truth[:, 0] = markerData[:, i]
        for k in range(0, 3, 1):
            pVector = np.zeros([tissueNum, 1])
            pVector[:, 0] = p0[k, :, i]
            temp = np.dot(reference, pVector)
            error = truth - temp
            error = np.sqrt(np.mean(np.multiply(error, error)))
            if error < minRMSE:
                minRMSE = error
                proportionDeconvolution[:, i] = p0[k, :, i]
                bestNu = k
        nuHistory.append(nu[bestNu])
    return proportionDeconvolution


def pre_marker_select(reference, markerData):
    """
    marker selection according to 1 vs 1 rule
    :param reference:
    :param markerData:
    :return:
    """
    cellTypeNumber = np.size(reference, 1)
    markerNumber = np.size(reference, 0)
    selectedCpG = np.arange(markerNumber)
    area = np.zeros(markerNumber)
    area = area.astype(np.bool)
    for i in range(cellTypeNumber):
        for j in range(i + 1, cellTypeNumber):
            temp = reference[:, [i, j]]
            tempSum = np.sum(temp, axis=1)
            tempSum1 = np.zeros([markerNumber, 2])
            tempSum1[:, 0] = tempSum
            tempSum1[:, 1] = tempSum
            temp = temp / tempSum1
            pairSortIndexIncrease = np.argsort(temp[:, 0])
            pairSortIndexDecrease = np.argsort(temp[:, 1])
            area[pairSortIndexIncrease[0 : 100]] = True
            area[pairSortIndexDecrease[0 : 100]] = True
    selectedCpG = selectedCpG[area]
    ref = reference[selectedCpG, :]
    mix = markerData[selectedCpG, :]
    return ref, mix


def filt_zeros(ref, mix):
    """
    remove markers where values of reference are zeros
    :param ref:
    :param mix:
    :return:
    """
    ref1 = []
    mix1 = []
    for i in range(np.size(ref, 0)):
        if np.max(ref[i, :]) > 0:
            ref1.append(ref[i, :])
            mix1.append(mix[i, :])
    return np.asarray(ref1), np.asarray(mix1)