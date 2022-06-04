# test pipeline for the bruyere dataset

from fn_cfg import *
import params as cfg

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
time = fileObjects[3]
trigChannel = fileObjects[4]

filtering = filters()
notchFilterOutput = filtering.notch(rawEEG,line,fs)
plots(time,notchFilterOutput,titles=cfg.channelNames,figsize=cfg.figure_size,pltclr=cfg.plot_color)

###############################################################################
original_EEG = notchFilterOutput[:,1]

"""
Detection of Ocular Artifacts using DWT (haar wavelet transform)
Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
and "ARDER: An Automatic EEG Artifacts Detection and Removal System"
"""
wavelet = 'haar'
coeffs = (wavedec(original_EEG,wavelet,level=5))

coeffs_approx = coeffs[0]
coeffs_5l = coeffs[1]
coeffs_4l = coeffs[2]
coeffs_3l = coeffs[3]
coeffs_2l = coeffs[4]
coeffs_1l = coeffs[5]

coeffs_approx = np.abs(coeffs_approx)
coeffs_5l = np.abs(coeffs_5l)
coeffs_4l = np.abs(coeffs_4l)
coeffs_3l = np.abs(coeffs_3l)
coeffs_2l = np.abs(coeffs_2l)
coeffs_1l = np.abs(coeffs_1l)

# Combine all extracted coefficients into the accepted coefficients array
acc_coeffs = []
acc_coeffs.append(coeffs_approx)
acc_coeffs.append(coeffs_5l)
acc_coeffs.append(coeffs_4l)
acc_coeffs.append(coeffs_3l)
acc_coeffs.append(coeffs_2l)
acc_coeffs.append(coeffs_1l)

# using the accepted coefficients, reconstruct the signal
# zones_OZ = areas with ocular artifacts
recov_signal = waverec(acc_coeffs,wavelet)
recov_signal[abs(recov_signal) > 100] = 0
zones_OZ = np.where(recov_signal == 0)[0]

# original signal with ocular artifacts
ocularArtifacts = original_EEG[zones_OZ]

plt.plot(recov_signal)
plt.title('clean signal')
plt.show()
plt.plot(ocularArtifacts)
plt.title('Ocular Artifacts')
plt.show()


"""
Removal of Ocular Artifacts
Based on "ARDER: An Automatic EEG Artifacts Detection and Removal System"
"""

wavelet = 'bior4.4'
coeffs = (wavedec(ocularArtifacts,wavelet,level=8))

# Extract approximate and detail coefficients
coeffs_approx = coeffs[0]
coeffs_8l = coeffs[1]
coeffs_7l = coeffs[2]
coeffs_6l = coeffs[3]
coeffs_5l = coeffs[4]
coeffs_4l = coeffs[5]
coeffs_3l = coeffs[6]
coeffs_2l = coeffs[7]
coeffs_1l = coeffs[8]

# Wavelet Thresholding for Denoising
# Based on Statistical Threshold (ST) method from 
# "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
threshold_approx = 1.5*np.std(coeffs_approx)
threshold_coeff8L = 1.5*(np.std(coeffs_8l))
threshold_coeff7L = 1.5*(np.std(coeffs_7l))
threshold_coeff6L = 1.5*(np.std(coeffs_6l))
threshold_coeff5L = 1.5*(np.std(coeffs_5l))
threshold_coeff4L = 1.5*(np.std(coeffs_4l))
threshold_coeff3L = 1.5*(np.std(coeffs_3l))
threshold_coeff2L = 1.5*(np.std(coeffs_2l))
threshold_coeff1L = 1.5*(np.std(coeffs_1l))

coeffs_approx[abs(coeffs_approx) > threshold_coeff8L] = 0
coeffs_8l[abs(coeffs_8l) > threshold_coeff8L] = 0
coeffs_7l[abs(coeffs_7l) > threshold_coeff8L] = 0
coeffs_6l[abs(coeffs_6l) > threshold_coeff8L] = 0
coeffs_5l[abs(coeffs_5l) > threshold_coeff8L] = 0
coeffs_4l[abs(coeffs_4l) > threshold_coeff8L] = 0
coeffs_3l[abs(coeffs_3l) > threshold_coeff8L] = 0
coeffs_2l[abs(coeffs_2l) > threshold_coeff8L] = 0
coeffs_1l[abs(coeffs_1l) > threshold_coeff8L] = 0

acc_coeffs = []
acc_coeffs.append(coeffs_approx)
acc_coeffs.append(coeffs_8l)
acc_coeffs.append(coeffs_7l)
acc_coeffs.append(coeffs_6l)
acc_coeffs.append(coeffs_5l)
acc_coeffs.append(coeffs_4l)
acc_coeffs.append(coeffs_3l)
acc_coeffs.append(coeffs_2l)
acc_coeffs.append(coeffs_1l)

# using the accepted coefficients, reconstruct the signal
# zones_OZ = areas with ocular artifacts
denoised_ocularArtifacts = waverec(acc_coeffs,wavelet)

plt.plot(denoised_ocularArtifacts,color='red')
plt.title('processed ocular artifacts')
plt.show()
plt.plot(ocularArtifacts,color='blue')
plt.show()

# Signal to Artifact Ratio (SAR)
#sar = 10*np.log10(np.std(ocularArtifacts)/np.std(ocularArtifacts-denoised_ocularArtifacts))
#print('Signal to Artifact Ratio (SAR):',sar)

# Mean Square Error (MSE)
# Smaller values of MSE means that the denoised EEG is more closed to clean EEG
#mse = np.mean((denoised_ocularArtifacts-ocularArtifacts)**2)/len(ocularArtifacts)
#print('Mean Squared Error (MSE):',mse)

# Dump denoised OA signal inside original EEG
s = original_EEG
l = zones_OZ
m = denoised_ocularArtifacts
d = dict(zip(l, m))
oa_cleanEEG = [d.get(i, j) for i, j in enumerate(s)]

# Signal to Artifact Ratio (SAR)
sar = 10*np.log10(np.std(original_EEG)/np.std(original_EEG-oa_cleanEEG))
print('Signal to Artifact Ratio (SAR):',sar)

# Mean Square Error (MSE)
# Smaller values of MSE means that the denoised EEG is more closed to clean EEG
mse = np.mean((oa_cleanEEG-original_EEG)**2)/len(original_EEG)
print('Mean Squared Error (MSE):',mse)

plt.plot(original_EEG,color='red',label='original')
plt.title('Original EEG')
plt.show()
plt.plot(oa_cleanEEG,color='blue',label='clean')
plt.title('Clean EEG')
plt.show()


def detect_remove_MA(original_EEG,fs):
    """
    -   Detection of Muscular Artifacts 
        -   Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
            and "ARDER: An Automatic EEG Artifacts Detection and Removal System"
        -   MA is detected if mad is greater than or equal to 1, hence removal tecnique is applied
    -   Removal of Muscular Artifacts 
        Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
        and "ARDER: An Automatic EEG Artifacts Detection and Removal System"
    """

    #   Detection of Muscular Artifacts 
    #   Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
    #   and "ARDER: An Automatic EEG Artifacts Detection and Removal System"
    #   MA is detected if mad is greater than or equal to 1, hence removal tecnique is applied
    freqs,psd = signal.welch(original_EEG,fs,nperseg=1024)

    # Extract Power Spectral Density (PSD) for 30-100 Hz and 0-100 Hz
    def find_nearest(X, value):
        return X[np.unravel_index(np.argmin(np.abs(X - value)), X.shape)]

    nearest_0Hz = find_nearest(freqs,0)
    nearest_30Hz = find_nearest(freqs,30)
    nearest_100Hz = find_nearest(freqs,100)

    idx_0Hz = np.where(freqs==nearest_0Hz)
    idx_30Hz = np.where(freqs==nearest_30Hz)
    idx_100Hz = np.where(freqs==nearest_100Hz)

    psd_30_100Hz = psd[idx_30Hz[0][0]:idx_100Hz[0][0]]
    psd_0_100Hz = psd[idx_0Hz[0][0]:idx_100Hz[0][0]]

    mad = np.sum(psd_30_100Hz)/np.sum(psd_0_100Hz)
    print('MAD Value:',mad)

    if mad > 0:

        #   Removal of Muscular Artifacts 
        #   Based on "Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data"
        #   and "ARDER: An Automatic EEG Artifacts Detection and Removal System"

        # Extract imfs from original EEG
        imfs = emd.sift.sift(original_EEG)
        x_t = imfs
        y_t = x_t-1
        cca = CCA(n_components=len(x_t.T))
        cca.fit(x_t,y_t)
        x_c,y_c= cca.transform(x_t,y_t)

        # Apply threshold to the CCA coefficients
        x_c[(x_c) > mad] = 0
        # Inverse CCA to get the denoised signal
        x_c_old = cca.inverse_transform(x_c)
        # Take mean of the IMFs
        recov_signal = np.mean(x_c_old,axis=1)

        plt.plot(recov_signal,color='red')
        plt.title('MA Recovered signal')
        plt.show()
        plt.plot(original_EEG,color='blue')
        plt.title('Original signal')

    elif mad <= 0:
        print('No MA detected')
        recov_signal = original_EEG
        print('No muscular artifacts detected')
        recov_signal = original_EEG
    return recov_signal

test = detect_remove_MA(original_EEG,fs)