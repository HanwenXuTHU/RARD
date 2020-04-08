# -*- coding: utf-8 -*-
"""
Created on Wen Apr 8 10:18 2020

@author: Hanwen Xu

E-mail: hw-xu16@mails.tsinghua.edu.cn

New update: process data with unknown content
"""
import numpy as np
from tqdm import tqdm
from scipy import stats, optimize
from sklearn.svm import NuSVR

def DeconCCN(reference, markerData, unknown = False):
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
        t = iterativeWeightedSVR(bestReference, bestMarker, unknown=unknown)
        bestNumberHistory.append(bestNumber)
        proportionDeconvolution[:, i] = t[:, 0]
    return proportionDeconvolution


def mixedConditionNumber(referenceSelect, markerSelect):
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


def iterativeWeightedSVR(reference, markerData, unknown=False):
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(markerData, 1)
    markerNumber = np.size(reference, 0)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples]) if unknown == False else np.zeros([tissueNum + 1, numOfSamples])
    for i in tqdm(range(numOfSamples)):
        iterNumber = 10
        iterReference = reference
        iterMarkerData = markerData[:, i]
        # Threshold Searching
        mixture = optimize.nnls(iterReference, iterMarkerData)
        test = mixture[0]
        t = test / np.sum(test) if unknown == False else test
        c1 = np.zeros([np.size(iterReference, 0), 1])
        c1[:, 0] = iterMarkerData[:]
        t1 = np.zeros([tissueNum, 1])
        t1[:, 0] = t
        c2 = np.dot(iterReference, t1)
        error = np.abs(c2 - c1)
        error = error[:, 0]
        error = np.sort(error)
        index = int(np.size(error) * 0.8)
        threshold = error[index]
        for j in range(iterNumber):
            mixture = optimize.nnls(iterReference, iterMarkerData)
            test = mixture[0]
            t = test / np.sum(test) if unknown == False else test
            c1 = np.zeros([np.size(iterReference, 0), 1])
            c1[:, 0] = iterMarkerData[:]
            t1 = np.zeros([tissueNum, 1])
            t1[:, 0] = t
            c2 = np.dot(iterReference, t1)
            error = np.abs(c2 - c1)
            error = error[:, 0]
            iterSelected = np.arange(np.size(iterReference, 0))
            area = error < threshold
            iterSelected = iterSelected[area]
            iterReference = iterReference[iterSelected, :]
            iterMarkerData = iterMarkerData[iterSelected]
        singleMarker = np.zeros([np.size(iterReference, 0), 1])
        singleMarker[:, 0] = iterMarkerData
        t = nuSVR(iterReference, singleMarker, unknown=unknown)
        t = t[:, 0]
        if unknown == False:
            proportionDeconvolution[:, i] = t
        else:
            proportionDeconvolution[0 : -1, i] = t
            proportionDeconvolution[-1, i] = max(1 - np.sum(t), 0)
    return proportionDeconvolution


def nuSVR(reference, markerData, unknown=False):
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
            tc = np.dot(clf.dual_coef_, clf.support_vectors_)
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
