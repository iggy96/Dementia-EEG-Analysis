# test pipeline for the bruyere dataset

from fn_cfg import *

localPath = '/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset'
filename = '0002_2_12122019_1225'

version = 1.0
#filename = '0_1_12072018_1206'
#localPath = '/Users/joshuaighalo/Downloads/raw'
y_lim = 1000
figsize = (8,8)
plot_color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
dispIMG = True
titles = ['Fz','Cz','Pz','P07','OZ']

device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]
rawEOG = fileObjects[1]
rawEEGEOG = fileObjects[2]
time_period = fileObjects[3]
trigChannel = fileObjects[4]

plots(time_period[0:50],rawEEGEOG[0:50],titles,figsize,plot_color)

import pywt
from pywt import wavedec, waverec
import numpy as np

x = rawEEG[:,1]
"""
Detection of Ocular Artifacts using DWT (haar wavelet transform)
Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
and "ARDER: An Automatic EEG Artifacts Detection and Removal System"
"""
wavelet = 'haar'
coeffs = (wavedec(x,wavelet,level=5))

coeffs_approx = coeffs[0]
coeffs_5l = coeffs[1]
coeffs_4l = coeffs[2]
coeffs_3l = coeffs[3]
coeffs_2l = coeffs[4]
coeffs_1l = coeffs[5]

# Wavelet Thresholding for Denoising
# Based on Statistical Threshold (ST) method from 
# "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
threshold_approx = 1.5*np.std(coeffs_approx)
threshold_coeff5L = 1.5*(np.std(coeffs_5l))
threshold_coeff4L = 1.5*(np.std(coeffs_4l))
threshold_coeff3L = 1.5*(np.std(coeffs_3l))
threshold_coeff2L = 1.5*(np.std(coeffs_2l))
threshold_coeff1L = 1.5*(np.std(coeffs_1l))

"""
use this when removing the artifacts
# Applying level 5 threshold
coeffs_approx[abs(coeffs_approx) > threshold_coeff5L] = 0
coeffs_5l[abs(coeffs_5l) > threshold_coeff5L] = 0
coeffs_4l[abs(coeffs_4l) > threshold_coeff5L] = 0
coeffs_3l[abs(coeffs_3l) > threshold_coeff5L] = 0
coeffs_2l[abs(coeffs_2l) > threshold_coeff5L] = 0
coeffs_1l[abs(coeffs_1l) > threshold_coeff5L] = 0
"""

coeffs_approx[abs(coeffs_approx) > threshold_approx] = 0
coeffs_5l[abs(coeffs_5l) > threshold_coeff5L] = 0
coeffs_4l[abs(coeffs_4l) > threshold_coeff4L] = 0
coeffs_3l[abs(coeffs_3l) > threshold_coeff3L] = 0
coeffs_2l[abs(coeffs_2l) > threshold_coeff2L] = 0
coeffs_1l[abs(coeffs_1l) > threshold_coeff1L] = 0

# Combine all extracted coefficients into the accepted coefficients array
acc_coeffs = []
acc_coeffs.append(coeffs_approx)
acc_coeffs.append(coeffs_5l)
acc_coeffs.append(coeffs_4l)
acc_coeffs.append(coeffs_3l)
acc_coeffs.append(coeffs_2l)
acc_coeffs.append(coeffs_1l)

###############################################################################
# using the accepted coefficients, reconstruct the signal
recov_signal = waverec(acc_coeffs,wavelet)
# zeros in the recov signal represent ocular artifacts
# mark parts of signal with amplitude > 100 as ocular artifacts
ocular_artifacts = np.where(abs(recov_signal) > 100)

plt.plot(time_period,recov_signal,color='red')
plt.show()
plt.plot(time_period,x,color='red')
plt.show()



