# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:00:13 2021

@author: oseho
"""


from df_lib import*
from params import*
import params as cfg


# Preprocessing functions
# Event Related Potentials (ERP) extracting functions 
class importFile:
    """
    Functionality is geared towards importing raw eeg files from the gtec system.
    The class contains a subclass which is the nuerocatch class
        -   The neurocatch class facilitates the selection of the version (either V1.0 or v1.1) of the neurocatch system
            utilized during the recording of the eeg data.
        -   This function also contains an EOG channel selector block which facilitates the correct selection of the 
            EOG channel.
        -   Function returns the raw EEG,raw EOG,an array holding both raw EEG and raw EOG,time (Ts) and the trigger
            channel (trig)
    """
    class neurocatch:
        def init(self,version,filename,localPath,dispIMG):
            if version == 1.0:
                data = filename
                localPath = localPath.replace(os.sep, '/')   
                localPath = localPath + '/'  
                path = localPath+data
                os.chdir(path)
                for file in (os.listdir(path)):
                    if file.endswith(".txt"):
                        file_path = f"{path}/{file}"
                        
                with open(file_path, 'r') as f:
                    lines = f.readlines(200)

                def two_digits(data):
                    # electrodes with Nan kOhms 
                    reject_impedance = 1000
                    fullstring = data
                    substring = "NaN"
                    if substring in fullstring:
                        val3 = reject_impedance
                    # electrodes with numeric kOhms
                    else:
                        val1 = data
                        val2 = val1[4:] # delete unnecessary characters from the front
                        val2 = val2[:2]
                        val3 = int(float(val2)) # delete unnecessary characters from the back then convert to integer
                        import math
                        # check 1
                        digits = int(math.log10(val3))+1 # check number of digits in result
                        if digits == 2: # expected result
                            val3 = val3
                        if digits < 2: # unexpected result 1
                            val5 = val1
                            val5 = val5[4:]
                            val5 = val5[:3]
                            val3 = int(float(val5))
                        # check 2
                        digits = int(math.log10(val3))+1 # check number of digits in result
                        if digits == 2: # expected result
                            val3 = val3
                        if digits < 1: # unexpected result 1
                            val6 = val1
                            val6 = val6[4:]
                            val6 = val6[:4]
                            val3 = int(float(val6))
                    return val3
                    
                def three_digits(data):
                    # electrodes with Nan kOhms 
                    reject_impedance = 1000
                    fullstring = data
                    substring = "NaN"
                    if substring in fullstring:
                        val3 = reject_impedance
                    # electrodes with numeric kOhms
                    else:
                        val1 = data
                        val2 = val1[4:]
                        val2 = val2[:3]
                        val3 = int(float(val2))
                        import math
                        # check 1
                        digits = int(math.log10(val3))+1 # check number of digits in result
                        if digits == 3: # expected result
                            val3 = val3
                        if digits < 3: # unexpected result 1
                            val5 = val1
                            val5 = val5[4:]
                            val5 = val5[:4]
                            val3 = int(float(val5))
                        # check 2
                        digits = int(math.log10(val3))+1 # check number of digits in result
                        if digits == 3: # expected result
                            val3 = val3
                        if digits < 2: # unexpected result 1
                            val6 = val1
                            val6 = val6[4:]
                            val6 = val6[:5]
                            val3 = int(float(val6))
                    return val3
                    
                # extract from metadata file, channels that collected eeg data
                device_chans = ['FZ','CZ','P3','PZ','P4','PO7','PO8','OZ','unknown channel','unknown channel','unknown channel','sampleNumber','battery','trigger']
                def lcontains(needle_s, haystack_l):
                    try: return [i for i in haystack_l if needle_s in i][0]
                    except IndexError: return None
                metadata_chans = [lcontains(device_chans[i],lines) for i in range(len(device_chans))]

                p3 = metadata_chans[2]
                if len(p3)<=15:
                    p3 = two_digits(p3)
                else:
                    p3 = three_digits(p3)
                p4 = metadata_chans[4]
                if len(p4)<=15:
                    p4 = two_digits(p4)
                else:
                    p4 = three_digits(p4)

                p07 = metadata_chans[5]
                if len(p07)<=15:
                    p07 = two_digits(p07)
                else:
                    p07 = three_digits(p07)

                p08 = metadata_chans[6]
                if len(p08)<=15:
                    p08 = two_digits(p08)
                else:
                    p08 = three_digits(p08)

                oz = metadata_chans[7]
                if len(oz)<=15:
                    oz = two_digits(oz)
                else:
                    oz = three_digits(oz)

            elif version == 1.1:
                localPath = localPath.replace(os.sep, '/')   
                localPath = localPath + '/'  
                path = localPath+filename
                os.chdir(path)
                for file in os.listdir():
                    if file.endswith(".json"):
                        file_path = f"{path}/{file}"

                metadata = open(file_path)
                metadata = json.load(metadata)
                if dispIMG == True:
                    for i in metadata:
                        print(i)
                    print(metadata['version'])
                else:
                    pass
                metadata_chans = metadata['channels']
                metadata_imp = metadata['impedances']
                # p3,p4,p07,p08,oz
                p3 = (metadata_imp[0])['P3']
                p4 = (metadata_imp[0])['P4']
                p07 = (metadata_imp[0])['PO7']
                p08 = (metadata_imp[0])['PO8']
                oz = (metadata_imp[0])['OZ']

            #  import raw file 1
            pathBin = [path+'/'+filename+'.bin']
            filenames = glob.glob(pathBin[0])
            if dispIMG == True:
                print(filenames)
            else:
                pass
            data = [np.fromfile(f, dtype=np.float32) for f in filenames]
            data1 = data[0]
            dataCols = len(metadata_chans)
            dataRows = int(len(data1)/dataCols)           
            data1 = data1.reshape(dataRows, dataCols)

            #  configure eeg channels
            eegChans = gtec['eegChans']
            fz = eegChans[0]
            fz1 = data1[:,fz]     #extract column 0
            fz = fz1.reshape(1,dataRows)

            cz = eegChans[1]
            cz1 = data1[:,cz]
            cz = cz1.reshape(1, dataRows)

            pz = eegChans[2]
            pz1 = data1[:,pz]
            pz = pz1.reshape(1,dataRows)

            # %% configure eog channels
            if p3 < 501:
                eogChans = gtec['eogChans']
                eog10 = eogChans[0]
                eogChan1 = data1[:,eog10]
                eogChan1 = eogChan1.reshape(1,dataRows)
                if dispIMG == True:
                    print('channel P3 utilized')
                else:
                    pass
            else:
                eogChan1 = np.zeros(len(fz.T))
                eogChan1 = eogChan1.reshape(1,len(eogChan1))

            if p4 < 501:
                eogChans = gtec['eogChans']
                eog20 = eogChans[1]
                eogChan2 = data1[:,eog20]
                eogChan2 = eogChan2.reshape(1,dataRows)
                if dispIMG == True:
                    print('channel P4 utilized')
                else:
                    pass
            else:
                eogChan2 = np.zeros(len(fz.T))
                eogChan2 = eogChan2.reshape(1,len(eogChan2))

            if p07 < 501:
                eogChans = gtec['eogChans']
                eog30 = eogChans[2]
                eogChan3 = data1[:,eog30]
                eogChan3 = eogChan3.reshape(1,dataRows)
                if dispIMG == True:
                    print('channel P07 utilized')
                else:
                    pass
            else:
                eogChan3 = np.zeros(len(fz.T))
                eogChan3 = eogChan3.reshape(1,len(eogChan3))

            if p08 < 501:
                eogChans = gtec['eogChans']
                eog40 = eogChans[3]
                eogChan4 = data1[:,eog40]
                eogChan4 = eogChan4.reshape(1,dataRows)
                if dispIMG == True:
                    print('channel P08 utilized')
                else:
                    pass
            else:
                eogChan4 = np.zeros(len(fz.T))
                eogChan4 = eogChan4.reshape(1,len(eogChan4))

            if oz < 501:
                eogChans = gtec['eogChans']
                eog50 = eogChans[4]
                eogChan5 = data1[:,eog50]
                eogChan5 = eogChan5.reshape(1,dataRows)
                if dispIMG == True:
                    print('channel 0Z utilized')
                else:
                    pass
            else:
                eogChan5 = np.zeros(len(fz.T))
                eogChan5 = eogChan5.reshape(1,len(eogChan5))

            # %% configure trigger channel
            trigCol = gtec['trigCol']
            trig = data1[:,trigCol]
            trig = trig.reshape(1,dataRows)

            # %% configure raw file
            rawData = np.concatenate((fz, cz, pz,eogChan1,eogChan2,eogChan3,eogChan4,eogChan5,trig))
            rawData = rawData.T

            # delete non zero columns i.e., the eogchans that are not in the data represented by zero columns
            mask = (rawData == 0).all(0)
                # Find the indices of these columns
            column_indices = np.where(mask)[0]
                # Update x to only include the columns where non-zero values occur.
            rawData = rawData[:,~mask] # rawData containing eegChans,

            # %% the new raw data just containing the required eegchans,eogchans and the trig channel
                # correctly name channels in raw Data
            csm = dict(fz=[0], cz=[1], pz=[2], eog1=[3], eog2=[4], eog3=[5], eog4=[6],eog5=[7], ntrig=[8])
            csm_fz = csm['fz']
            csm_fz = csm_fz[0]
            csm_cz = csm['cz']
            csm_cz = csm_cz[0]
            csm_pz = csm['pz']
            csm_pz = csm_pz[0]
            csm_eog1 = csm['eog1']
            csm_eog1 = csm_eog1[0]
            csm_eog2 = csm['eog2']
            csm_eog2 = csm_eog2[0]
            csm_eog3 = csm['eog3']
            csm_eog3 = csm_eog3[0]
            csm_eog4 = csm['eog4']
            csm_eog4 = csm_eog4[0]
            csm_eog5 = csm['eog5']
            csm_eog5 = csm_eog5[0]
            csm_ntrig = csm['ntrig']
            csm_ntrig = csm_ntrig[0]

            if len(rawData.T)==4:
                csm_ntrig = 3
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_ntrig]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & no EOG channels')
                else:
                    pass
            
            elif len(rawData.T)==5:
                csm_ntrig = 4
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & one EOG channel')
                else:
                    pass

            elif len(rawData.T)==6:
                csm_ntrig = 5
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & two EOG channels')
                else:
                    pass

            elif len(rawData.T)==7:
                csm_ntrig = 6
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & three EOG channels')
                else:
                    pass

            elif len(rawData.T)==8:
                csm_ntrig = 7
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3,csm_eog4]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & four EOG channels')
                else:
                    pass

            elif len(rawData.T)==9:
                csm_ntrig = 8
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_eog5,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_eog5]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                if dispIMG == True:
                    print('data contains Fz, Cz, Pz & five EOG channels')
                else:
                    pass

            # time period of scan
            fs = gtec['fs']
            dt = 1/fs
            stop = dataRows/fs
            Ts = (np.arange(0,stop,dt)).reshape(len(np.arange(0,stop,dt)),1)
            return rawEEG,rawEOG,rawEEGEOG,Ts,trig

class filters:
    """
     filters for EEG data
     filtering order: adaptive filter -> notch filter -> bandpass filter (or lowpass filter, highpass filter)
    """
    def notch(self,data,line,fs,Q=30):
        """
           Inputs  :   data    - 2D numpy array (d0 = samples, d1 = channels) of unfiltered EEG data
                       cut     - frequency to be notched (defaults to config)
                       fs      - sampling rate of hardware (defaults to config)
                       Q       - Quality Factor (defaults to 30) that characterizes notch filter -3 dB bandwidth bw relative to its center frequency, Q = w0/bw.   
           Output  :   y     - 2D numpy array (d0 = samples, d1 = channels) of notch-filtered EEG data
           NOTES   :   
           Todo    : report testing filter characteristics
        """
        cut = line
        w0 = cut/(fs/2)
        b, a = signal.iirnotch(w0, Q)
        y = signal.filtfilt(b, a, data, axis=0)
        return y

    def butterBandPass(self,data,lowcut,highcut,fs,order=4):
        """
           Inputs  :   data    - 2D numpy array (d0 = samples, d1 = channels) of unfiltered EEG data
                       low     - lower limit in Hz for the bandpass filter (defaults to config)
                       high    - upper limit in Hz for the bandpass filter (defaults to config)
                       fs      - sampling rate of hardware (defaults to config)
                       order   - the order of the filter (defaults to 4)  
           Output  :   y     - 2D numpy array (d0 = samples, d1 = channels) of notch-filtered EEG data
           NOTES   :   
           Todo    : report testing filter characteristics
         data: eeg data (samples, channels)
         some channels might be eog channels
        """
        low_n = lowcut
        high_n = highcut
        sos = butter(order, [low_n, high_n], btype="bandpass", analog=False, output="sos",fs=fs)
        y = sosfiltfilt(sos, data, axis=0)
        return y

    def adaptive(self,eegData,eogData,nKernel=5, forgetF=0.995,  startSample=0, p = False):
        """
           Inputs:
           eegData - A matrix containing the EEG data to be filtered here each channel is a column in the matrix, and time
           starts at the top row of the matrix. i.e. size(data) = [numSamples,numChannels]
           eogData - A matrix containing the EOG data to be used in the adaptive filter
           startSample - the number of samples to skip for the calculation (i.e. to avoid the transient)
           p - plot AF response (default false)
           nKernel = Dimension of the kernel for the adaptive filter
           Outputs:
           cleanData - A matrix of the same size as "eegdata", now containing EOG-corrected EEG data.
           Adapted from He, Ping, G. Wilson, and C. Russell. "Removal of ocular artifacts from electro-encephalogram by adaptive filtering." Medical and biological engineering and computing 42.3 (2004): 407-412.
        """
        #   reshape eog array if necessary
        if len(eogData.shape) == 1:
            eogData = np.reshape(eogData, (eogData.shape[0], 1))
        # initialise Recursive Least Squares (RLS) filter state
        nEOG = eogData.shape[1]
        nEEG = eegData.shape[1]
        hist = np.zeros((nEOG, nKernel))
        R_n = np.identity(nEOG * nKernel) / 0.01
        H_n = np.zeros((nEOG * nKernel, nEEG))
        X = np.hstack((eegData, eogData)).T          # sort EEG and EOG channels, then transpose into row variables
        eegIndex = np.arange(nEEG)                              # index of EEG channels within X
        eogIndex = np.arange(nEOG) + eegIndex[-1] + 1           # index of EOG channels within X
        for n in range(startSample, X.shape[1]):
            hist = np.hstack((hist[:, 1:], X[eogIndex, n].reshape((nEOG, 1))))  # update the EOG history by feeding in a new sample
            tmp = hist.T                                                        # make it a column variable again (?)
            r_n = np.vstack(np.hsplit(tmp, tmp.shape[-1]))
            K_n = np.dot(R_n, r_n) / (forgetF + np.dot(np.dot(r_n.T, R_n), r_n))                                           # Eq. 25
            R_n = np.dot(np.power(forgetF, -1),R_n) - np.dot(np.dot(np.dot(np.power(forgetF, -1), K_n), r_n.T), R_n)       #Update R_n
            s_n = X[eegIndex, n].reshape((nEEG, 1))                   #get EEG signal and make sure it's a 1D column array
            e_nn = s_n - np.dot(r_n.T, H_n).T  #Eq. 27
            H_n = H_n + np.dot(K_n, e_nn.T)
            e_n = s_n - np.dot(r_n.T, H_n).T
            X[eegIndex, n] = np.squeeze(e_n)
        cleanData = X[eegIndex, :].T
        return cleanData

    def butter_lowpass(self,data,cutoff,fs,order):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
        y = signal.lfilter(b, a, data)
        return y

    def butter_highpass(self,data,cutoff,fs,order):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
        y = signal.filtfilt(b, a, data)
        return y

def rising_edge(data):
    # used in trigger channel development before epoching
    trig = data
    trg = (trig >= 1) & (trig < 9)
    trig[trg > 1] = 0     # cut off any onset offset triggers outside the normal range.
    diff = np.diff(trig)
    t = diff > 0  # rising edges assume true
    k = np.array([ False])
    k = k.reshape(1,len(k))
    pos = np.concatenate((k, t),axis=1)
    trig[~pos] = 0
    newtrig = trig 
    trigCol = newtrig
    trigCol = trigCol.T
    trigCol = trigCol.reshape(len(trigCol))
    return trigCol

def peaktopeak(data):
    # used for artifact rejection
    a = np.amax(data, axis = 1)
    b = np.amin(data, axis = 1)
    p2p = a-b
    return p2p

class erpExtraction:
    """
      Inputs: trigger data produced from rising_edge()
            : standard and deviant tones elicit the N100 and P300 erps
            : congruent and incongruent words elicit the N400 erps
      Outputs: ERP data (samples, channels) for N100,P300 and N400
      Notes:   - The trigger channel is assumed to be the last channel in the data matrix
    """
    def N100P300(self,chanNames,scanID,trigger_channel,eegData,period,stimTrig,clip,dispIMG):
        """
          Inputs: trigger channels, bandpass filtered data for all channels, time period,
                  stimTrig, clip value
        """
        chanNames,scanID,trigger_channel,channel_data,period,stimTrig,clip,dispIMG = chanNames,scanID,trigger_channel,eegData,period,stimTrig,clip,dispIMG
        def algorithm(chanNames,scanID,trigger_channel,channel_data,period,stimTrig,clip,dispIMG):
            trigger_data = rising_edge(trigger_channel)
            trigCol = trigger_data
            avg = channel_data 
            Ts = period
            # STANDARD TONE [1]: extract the time points where stimuli 1 exists
            std = stimTrig['std']
            std = std[0]
            
            no_ones = np.count_nonzero(trigCol==1)
            if dispIMG == True:
                print("number of std tone event codes:",no_ones)
            else:
                pass
            
            result = np.where(trigCol == std)
            idx_10 = result[0]
            idx_11 = idx_10.reshape(1,len(idx_10))
            
            # sort target values
            desiredCols = trigCol[idx_10]
            desiredCols = desiredCols.reshape(len(desiredCols),1)
            
            # sort values before target
            idx_st = (idx_11 - 50)
            i = len(idx_st.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = trigCol[int(idx_st[:,x]):int(idx_11[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            startVals = allArrays
            startCols = startVals.reshape(no_ones,50)
            
            # sort values immediately after target
            idx_en11 = idx_11 + 1
            postTargetCols = trigCol[idx_en11]
            postTargetCols = postTargetCols.T
            
            # sort end values
                # ----------------------------------------------------------------------------
                # determination of the number steps to reach the last point of each epoch
                # the event codes are not evenly timed hence the need for steps determination
            a = trigCol[idx_10]
            a = a.reshape(len(a),1)
            
            b = Ts[idx_10]
            
            c_i = float("%0.3f" % (Ts[int(len(Ts)-1)]))
            c = (np.array([0,c_i]))
            c = c.reshape(1,len(c))
            
            std_distr = np.concatenate((a, b),axis=1)
            std_distr = np.vstack((std_distr,c))
            
            std_diff = np.diff(std_distr[:,1],axis=0)
            std_step = ((std_diff/0.002)-1)
            std_step = np.where(std_step > 447, 447, std_step)
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(std_step == std_step[0])
            if result:
                # use for equal epoch steps
                # sort end values
                idx_en12 = idx_en11 + std_step.T
                i = len(idx_en12.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = trigCol[int(idx_en11[:,x]):int(idx_en12[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                endCols = endVals.reshape(no_ones,447)
                
                # merge the different sections of the epoch to form the epochs
                trig_std = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            else:
                # use when we have unequal epoch steps
                    # --------------------------------------------
                    # apply the step to get the index of the last point of the epoch
                idx_en12 = idx_en11 + std_step.T 
                i = len(idx_en12.T)
                allArrays = []
                for x in range(i):
                    myArray = trigCol[int(idx_en11[:,x]):int(idx_en12[:,x])]
                    allArrays.append(myArray)
                endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                endVals = pd.DataFrame(endVals)
                endVals.fillna(endVals.mean(),inplace=True)
                endCols = endVals.values
                
                # merge the different sections of the epoch to form the epochs
                trig_std = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            # ----------------------------------------------------------------------------------------------------------------
            # implement trig epoch creation to eeg data
                # replace trigCol with avg 
                # get start, desired, post and end col, then merge together
                # remove data equal to non 1 or o triggers
            
            # sort target values
            eeg_desiredCols = avg[idx_10]
            eeg_desiredCols = eeg_desiredCols.reshape(len(eeg_desiredCols),1)
            
            # sort values before target
            i = len(idx_st.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = avg[int(idx_st[:,x]):int(idx_11[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            eeg_startVals = allArrays
            eeg_startCols = eeg_startVals.reshape(no_ones,50)
            
            # sort values immediately after target
            eeg_postTargetCols = avg[idx_en11]
            eeg_postTargetCols = eeg_postTargetCols.T
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(std_step == std_step[0])
            if result:
                # use for equal epochs
                # sort end values
                idx_en12 = idx_en11 + std_step.T
                i = len(idx_en12.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = avg[int(idx_en11[:,x]):int(idx_en12[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                eeg_endCols = endVals.reshape(no_ones,447)
                
                # merge the different sections of the epoch to form the epochs
                # eeg_std = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                epochs_std = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                
            else:
                # use for unequal epochs
                # sort end values
                idx_en12 = idx_en11 + std_step.T 
                i = len(idx_en12.T)
                allArrays = []
                for x in range(i):
                    myArray = avg[int(idx_en11[:,x]):int(idx_en12[:,x])]
                    allArrays.append(myArray)
                eeg_endVals = allArrays

                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                eeg_endVals = pd.DataFrame(eeg_endVals)
                eeg_endVals.fillna(eeg_endVals.mean(),inplace=True)
                eeg_endCols = eeg_endVals.values
                
                # merge the different sections of the epoch to form the epochs
                epochs_std = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
            
            # baseline correction
            prestim = epochs_std[:,0:49]
            mean_prestim = np.mean(prestim,axis=1)
            mean_prestim = mean_prestim.reshape(len(mean_prestim),1)
            bc_std = epochs_std - mean_prestim
            
            # artefact rejection
            p2p = peaktopeak(bc_std)
            result = np.where(p2p > clip)
            row = result[0]
            ar_std = np.delete(bc_std,(row),axis = 0)
            dif = ((len(bc_std)-len(ar_std))/len(bc_std))*100
            if len(ar_std) == len(bc_std):
                if dispIMG == True:
                    print("notice! epochs lost for std tone:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            elif len(ar_std) < len(bc_std):
                if dispIMG == True:
                    print("callback! epochs lost for std tone:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            
            if ar_std.size == 0:
                avg_std = np.zeros(len(bc_std.T))
                print(chanNames,':',scanID," removed for standard tones N1P3 analysis as all its epochs exceed the clip value of",clip)
            elif ar_std.size != 0:
                # averaging
                avg_std = np.mean(ar_std,axis=0)
            avg_std = avg_std

            #%%
            # DEVIANT TONE [2]: extract the time points where stimuli 2 exists
            dev = stimTrig['dev']
            dev = dev[0]
            
            no_twos = np.count_nonzero(trigCol==2)
            if dispIMG == True:
                print("number of dev tone event codes:",no_twos)
            else:
                pass
            
            result = np.where(trigCol == dev)
            idx_20 = result[0]
            idx_21 = idx_20.reshape(1,len(idx_20))
            
            # sort target values
            desiredCols = trigCol[idx_20]
            desiredCols = desiredCols.reshape(len(desiredCols),1)
            
            # sort values before target
            idx_dev = (idx_21 - 50)
            i = len(idx_dev.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = trigCol[int(idx_dev[:,x]):int(idx_21[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            startVals = allArrays
            startCols = startVals.reshape(no_twos,50)
            
            # sort values immediately after target
            idx_en21 = idx_21 + 1
            postTargetCols = trigCol[idx_en21]
            postTargetCols = postTargetCols.T
            
            # sort end values
                # ----------------------------------------------------------------------------
                # determination of the number steps to reach the last point of each epoch
                # the event codes are not evenly timed hence the need for steps determination
            a = trigCol[idx_20]
            a = a.reshape(len(a),1)
            
            b = Ts[idx_20]
            
            c_i = float("%0.3f" % (Ts[int(len(Ts)-1)]))
            c = (np.array([0,c_i]))
            c = c.reshape(1,len(c))
            
            dev_distr = np.concatenate((a, b),axis=1)
            dev_distr = np.vstack((dev_distr,c))
            
            dev_diff = np.diff(dev_distr[:,1],axis=0)
            dev_step = ((dev_diff/0.002)-1)
            dev_step = np.where(dev_step > 447, 447, dev_step)
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(dev_step == dev_step[0])
            if result:
                # sort end values
                idx_en22 = idx_en21 + dev_step.T
                i = len(idx_en22.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = trigCol[int(idx_en21[:,x]):int(idx_en22[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                endCols = endVals.reshape(no_twos,447)
            
                    # merge the different sections of the epoch to form the epochs
                trig_dev = np.concatenate((startCols,desiredCols,endCols),axis=1)
            else:
                # use when we have unequal epoch steps
                    # apply the step to get the index of the last point of the epoch
                idx_en22 = idx_en21 + dev_step.T 
                i = len(idx_en22.T)
                allArrays = []
                for x in range(i):
                    myArray = trigCol[int(idx_en21[:,x]):int(idx_en22[:,x])]
                    allArrays.append(myArray)
                endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                endVals = pd.DataFrame(endVals)
                endVals.fillna(endVals.mean(),inplace=True)
                endCols = endVals.values
                
                    # merge the different sections of the epoch to form the epochs
                trig_devpad = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            # ----------------------------------------------------------------------------------------------------------------
            # implement trig epoch creation to eeg data
                # replace trigCol with avg 
                # get start, desired, post and end col, then merge together
                # remove data equal to non 1 or o triggers
            
            # sort target values
            eeg_desiredCols = avg[idx_20]
            eeg_desiredCols = eeg_desiredCols.reshape(len(eeg_desiredCols),1)
            
            # sort values before target
            i = len(idx_dev.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = avg[int(idx_dev[:,x]):int(idx_21[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            eeg_startVals = allArrays
            eeg_startCols = eeg_startVals.reshape(no_twos,50)
            
            # sort values immediately after target
            eeg_postTargetCols = avg[idx_en21]
            eeg_postTargetCols = eeg_postTargetCols.T
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(dev_step == dev_step[0])
            if result:
                # use for equal epochs
                # sort end values
                idx_en22 = idx_en21 + dev_step.T
                i = len(idx_en22.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = avg[int(idx_en21[:,x]):int(idx_en22[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                eeg_endCols = endVals.reshape(no_twos,447)
                
                # merge the different sections of the epoch to form the epochs
                # eeg_dev = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                epochs_dev = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                
            else:
                # use for unequal epochs: fill steps with zeros 
                # sort end values
                idx_en22 = idx_en21 + dev_step.T 
                i = len(idx_en22.T)
                allArrays = []
                for x in range(i):
                    myArray = avg[int(idx_en21[:,x]):int(idx_en22[:,x])]
                    allArrays.append(myArray)
                eeg_endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                eeg_endVals = pd.DataFrame(eeg_endVals)
                eeg_endVals.fillna(eeg_endVals.mean(),inplace=True)
                eeg_endCols = eeg_endVals.values
                
                # merge the different sections of the epoch to form the epochs
                epochs_dev = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                # eeg_dev = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
            
            # baseline correction
            prestim = epochs_dev[:,0:49]
            mean_prestim = np.mean(prestim,axis=1)
            mean_prestim = mean_prestim.reshape(len(mean_prestim),1)
            bc_dev = epochs_dev - mean_prestim
            
            # artefact rejection
            p2p = peaktopeak(bc_dev)
            result = np.where(p2p > clip)
            row = result[0]
            ar_dev = np.delete(bc_dev,(row),axis = 0)
            dif = ((len(bc_dev)-len(ar_dev))/len(bc_dev))*100
            if len(ar_dev) == len(bc_dev):
                if dispIMG == True:
                    print("notice! epochs lost for dev tone:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            elif len(ar_dev) < len(bc_dev):
                if dispIMG == True:
                    print("callback! epochs lost for dev tone:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            
            if ar_dev.size == 0:
                avg_dev = np.zeros(len(bc_dev.T))
                print(chanNames,':',scanID," removed for deviant tones N1P3 analysis as all its epochs exceed the clip value of",clip)
            elif ar_dev.size != 0:
                avg_dev = np.mean(ar_dev,axis=0)
            avg_dev = avg_dev
            return avg_std,avg_dev,ar_std,ar_dev
        # algorithm function returns avg_std,avg_dev,ar_std,ar_dev
        out_final = []
        for i in range(len(channel_data.T)):
            out_final.append(algorithm(chanNames[i],scanID,trigger_channel,channel_data[:,i],period,stimTrig,clip,dispIMG))
        out_final = np.asarray(out_final,dtype="object").T
        out_final = out_final.transpose()
        return out_final

    def N400(self,chanNames,scanID,trigger_channel,eegData,period,stimTrig,clip,dispIMG):
        """
          Inputs: trigger channels, bandpass filtered data for all channels, time period,
                  stimTrig, clip value
        """
        
        chanNames,scanID,trigger_channel,channel_data,period,stimTrig,clip,dispIMG = chanNames,scanID,trigger_channel,eegData,period,stimTrig,clip,dispIMG
        def algorithm(chanNames,scanID,trigger_channel,channel_data,period,stimTrig,clip,dispIMG):
            trigger_data = rising_edge(trigger_channel)
            trigCol = trigger_data
            avg = channel_data 
            Ts = period
            # congruent word [4,7]: extract the time points where stimuli 1 exists
            con = stimTrig['con']
            con = con[0:2]
            con = (np.array([con])).T
            
            no_fours = np.count_nonzero(trigCol==4)
            no_sevens = np.count_nonzero(trigCol==7)
            no_cons = no_fours + no_sevens
            if dispIMG == True:
                print("number of con word event codes:",no_cons)
            else:
                pass
            
            result = np.where(trigCol == con[0])
            idx_30i = result[0]
            idx_30i = idx_30i.reshape(len(idx_30i),1)
            result = np.where(trigCol == con[1])
            idx_30ii = result[0]
            idx_30ii = idx_30ii.reshape(len(idx_30ii),1)
            idx_30 = np.vstack((idx_30i,idx_30ii))
            idx_31 = idx_30.reshape(1,len(idx_30))
            
            
            # sort target values
            desiredCols = trigCol[idx_30]
            desiredCols = desiredCols.reshape(len(desiredCols),1)
            
            # sort values before target
            idx_con = (idx_31 - 50)
            i = len(idx_con.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = trigCol[int(idx_con[:,x]):int(idx_31[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            startVals = allArrays
            startCols = startVals.reshape(no_cons,50)
            
            # sort values immediately after target
            idx_en31 = idx_31 + 1
            postTargetCols = trigCol[idx_en31]
            postTargetCols = postTargetCols.T
            
            # sort end values
                # ----------------------------------------------------------------------------
                # determination of the number steps to reach the last point of each epoch
                # the event codes are not evenly timed hence the need for steps determination
            a = trigCol[idx_30]
            a = a.reshape(len(a),1)
            
            b = Ts[idx_30]
            b = b.reshape(len(idx_30),1)
            
            c_i = float("%0.3f" % (Ts[int(len(Ts)-1)]))
            c = (np.array([0,c_i]))
            c = c.reshape(1,len(c))
            
            con_distr = np.concatenate((a, b),axis=1)
            con_distr = np.vstack((con_distr,c))
            
            con_diff = np.diff(con_distr[:,1],axis=0)
            con_step = ((con_diff/0.002)-1)
            con_step = np.where(con_step > 447, 447, con_step)
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(con_step == con_step[0])
            if result:
                # use for equal epoch steps
                # sort end values
                idx_en32 = idx_en31 + con_step.T
                i = len(idx_en32.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = trigCol[int(idx_en31[:,x]):int(idx_en32[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                endCols = endVals.reshape(no_cons,447)
                
                # merge the different sections of the epoch to form the epochs
                trig_con = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            else:
                # use when we have unequal epoch steps
                    # --------------------------------------------
                    # apply the step to get the index of the last point of the epoch
                idx_en32 = idx_en31 + con_step.T 
                i = len(idx_en32.T)
                allArrays = []
                for x in range(i):
                    myArray = trigCol[int(idx_en31[:,x]):int(idx_en32[:,x])]
                    allArrays.append(myArray)
                endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                endVals = pd.DataFrame(endVals)
                endVals.fillna(endVals.mean(),inplace=True)
                endCols = endVals.values
                
                # merge the different sections of the epoch to form the epochs
                trig_conpad = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            # ----------------------------------------------------------------------------------------------------------------
            # implement trig epoch creation to eeg data
                # replace trigCol with avg 
                # get start, desired, post and end col, then merge together
                # remove data equal to non 1 or o triggers
            
            # sort target values
            eeg_desiredCols = avg[idx_30]
            eeg_desiredCols = eeg_desiredCols.reshape(len(eeg_desiredCols),1)
            
            # sort values before target
            i = len(idx_con.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = avg[int(idx_con[:,x]):int(idx_31[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            eeg_startVals = allArrays
            eeg_startCols = eeg_startVals.reshape(no_cons,50)
            
            # sort values immediately after target
            eeg_postTargetCols = avg[idx_en31]
            eeg_postTargetCols = eeg_postTargetCols.T
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(con_step == con_step[0])
            if result:
                # use for equal epochs
                # sort end values
                idx_en32 = idx_en31 + con_step.T
                i = len(idx_en32.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = avg[int(idx_en31[:,x]):int(idx_en32[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                eeg_endCols = endVals.reshape(no_cons,447)
                
                # merge the different sections of the epoch to form the epochs
                # eeg_std = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                epochs_con = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                
            else:
                # use for unequal epochs
                # sort end values
                idx_en32 = idx_en31 + con_step.T 
                i = len(idx_en32.T)
                allArrays = []
                for x in range(i):
                    myArray = avg[int(idx_en31[:,x]):int(idx_en32[:,x])]
                    allArrays.append(myArray)
                eeg_endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                eeg_endVals = pd.DataFrame(eeg_endVals)
                eeg_endVals.fillna(eeg_endVals.mean(),inplace=True)
                eeg_endCols = eeg_endVals.values
                
                # merge the different sections of the epoch to form the epochs
                epochs_con = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
            
            
            # baseline correction
            prestim = epochs_con[:,0:49]
            mean_prestim = np.mean(prestim,axis=1)
            mean_prestim = mean_prestim.reshape(len(mean_prestim),1)
            bc_con = epochs_con - mean_prestim
            
            # artefact rejection
            p2p = peaktopeak(bc_con)
            result = np.where(p2p > clip)
            row = result[0]
            ar_con = np.delete(bc_con,(row),axis = 0)
            dif = ((len(bc_con)-len(ar_con))/len(bc_con))*100
            if len(ar_con) == len(bc_con):
                if dispIMG == True:
                    print("notice! epochs lost for con word:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            elif len(ar_con) < len(bc_con):
                if dispIMG == True:
                    print("callback! epochs lost for con word:","{:.2%}".format((int(dif))/100))
                else:
                    pass

            if ar_con.size == 0:
                avg_con = np.zeros(len(bc_con.T))
                print(chanNames,':',scanID," removed for congruent words N4 analysis as all its epochs exceed the clip value of",clip)
            else:
                # averaging
                avg_con = np.mean(ar_con,axis=0)
            avg_con = avg_con
            
            
            # %%
            # incongruent word [4,7]: extract the time points where stimuli 1 exists
            inc = stimTrig['inc']
            inc = inc[0:2]
            inc = (np.array([inc])).T
            
            no_fives = np.count_nonzero(trigCol==5)
            no_eights = np.count_nonzero(trigCol==8)
            no_incs = no_fives + no_eights
            if dispIMG == True:
                print("number of inc word event codes:",no_incs)
            else:
                pass
            
            result = np.where(trigCol == inc[0])
            idx_40i = result[0]
            idx_40i = idx_40i.reshape(len(idx_40i),1)
            result = np.where(trigCol == inc[1])
            idx_40ii = result[0]
            idx_40ii = idx_40ii.reshape(len(idx_40ii),1)
            idx_40 = np.vstack((idx_40i,idx_40ii))
            idx_41 = idx_40.reshape(1,len(idx_40))
            
            
            # sort target values
            desiredCols = trigCol[idx_40]
            desiredCols = desiredCols.reshape(len(desiredCols),1)
            
            # sort values before target
            idx_inc = (idx_41 - 50)
            i = len(idx_inc.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = trigCol[int(idx_inc[:,x]):int(idx_41[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            startVals = allArrays
            startCols = startVals.reshape(no_incs,50)
            
            # sort values immediately after target
            idx_en41 = idx_41 + 1
            postTargetCols = trigCol[idx_en41]
            postTargetCols = postTargetCols.T
            
            # sort end values
                # ----------------------------------------------------------------------------
                # determination of the number steps to reach the last point of each epoch
                # the event codes are not evenly timed hence the need for steps determination
            a = trigCol[idx_40]
            a = a.reshape(len(a),1)
            
            b = Ts[idx_40]
            b = b.reshape(len(idx_40),1)
            
            c_i = float("%0.3f" % (Ts[int(len(Ts)-1)]))
            c = (np.array([0,c_i]))
            c = c.reshape(1,len(c))
            
            inc_distr = np.concatenate((a, b),axis=1)
            inc_distr = np.vstack((inc_distr,c))
            
            inc_diff = np.diff(inc_distr[:,1],axis=0)
            inc_step = ((inc_diff/0.002)-1)
            inc_step = np.where(inc_step > 447, 447, inc_step)
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(inc_step == inc_step[0])
            if result:
                # use for equal epoch steps
                # sort end values
                idx_en42 = idx_en41 + inc_step.T
                i = len(idx_en42.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = trigCol[int(idx_en41[:,x]):int(idx_en42[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                endCols = endVals.reshape(no_inc,447)
                
                # merge the different sections of the epoch to form the epochs
                trig_inc = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            else:
                # use when we have unequal epoch steps
                    # --------------------------------------------
                    # apply the step to get the index of the last point of the epoch
                idx_en42 = idx_en41 + inc_step.T 
                i = len(idx_en42.T)
                allArrays = []
                for x in range(i):
                    myArray = trigCol[int(idx_en41[:,x]):int(idx_en42[:,x])]
                    allArrays.append(myArray)
                endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                l = endVals
                max_len = max([len(arr) for arr in l])
                padded = np.array([np.lib.pad(arr, (0, max_len - len(arr)), 'constant', constant_values=0) for arr in l])
                    # appropriate end columns
                endCols = padded
                
                # merge the different sections of the epoch to form the epochs
                trig_incpad = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
            
            # ----------------------------------------------------------------------------------------------------------------
            # implement trig epoch creation to eeg data
                # replace trigCol with avg 
                # get start, desired, post and end col, then merge together
                # remove data equal to non 1 or o triggers
            
            # sort target values
            eeg_desiredCols = avg[idx_40]
            eeg_desiredCols = eeg_desiredCols.reshape(len(eeg_desiredCols),1)
            
            # sort values before target
            i = len(idx_inc.T)
            allArrays = np.array([])
            for x in range(i):
                myArray = avg[int(idx_inc[:,x]):int(idx_41[:,x])]
                allArrays = np.concatenate([allArrays,myArray])
            eeg_startVals = allArrays
            eeg_startCols = eeg_startVals.reshape(no_incs,50)
            
            # sort values immediately after target
            eeg_postTargetCols = avg[idx_en41]
            eeg_postTargetCols = eeg_postTargetCols.T
            
            # check if the number of steps are the same if yes, meaning the epochs have the same length
            result = np.all(inc_step == inc_step[0])
            if result:
                # use for equal epochs
                # sort end values
                idx_en42 = idx_en41 + inc_step.T
                i = len(idx_en42.T)
                allArrays = np.array([])
                for x in range(i):
                    myArray = avg[int(idx_en41[:,x]):int(idx_en42[:,x])]
                    allArrays = np.concatenate([allArrays,myArray])
                endVals = allArrays
                eeg_endCols = endVals.reshape(no_incs,447)
                
                # merge the different sections of the epoch to form the epochs
                # eeg_std = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                epochs_inc = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
                
            else:
                # use for unequal epochs
                # sort end values
                idx_en42 = idx_en41 + inc_step.T 
                i = len(idx_en42.T)
                allArrays = []
                for x in range(i):
                    myArray = avg[int(idx_en41[:,x]):int(idx_en42[:,x])]
                    allArrays.append(myArray)
                eeg_endVals = allArrays
                    # add zeros to fill enable epochs of shorter length equal with epochs of adequate length
                eeg_endVals = pd.DataFrame(eeg_endVals)
                eeg_endVals.fillna(eeg_endVals.mean(),inplace=True)
                eeg_endCols = eeg_endVals.values
                
                # merge the different sections of the epoch to form the epochs
                epochs_inc = np.concatenate((eeg_startCols,eeg_desiredCols,eeg_endCols),axis=1)
            
            # baseline correction
            prestim = epochs_inc[:,0:49]
            mean_prestim = np.mean(prestim,axis=1)
            mean_prestim = mean_prestim.reshape(len(mean_prestim),1)
            bc_inc = epochs_inc - mean_prestim
            
            # artefact rejection
            p2p = peaktopeak(bc_inc)
            result = np.where(p2p > clip)
            row = result[0]
            ar_inc = np.delete(bc_inc,(row),axis = 0)
            dif = ((len(bc_inc)-len(ar_inc))/len(bc_inc))*100
            if len(ar_inc) == len(bc_inc):
                if dispIMG == True:
                    print("notice! epochs lost for inc word:","{:.2%}".format((int(dif))/100))
                else:
                    pass
            elif len(ar_inc) < len(bc_inc):
                if dispIMG == True:
                    print("callback! epochs lost for inc word:","{:.2%}".format((int(dif))/100))
                else:
                    pass

            if ar_inc.size == 0:
                avg_inc = np.zeros(len(bc_inc.T))
                print(chanNames,':',scanID," removed for incongruent words N4 analysis as all its epochs exceed the clip value of",clip)
            else:
                # averaging
                avg_inc = np.mean(ar_inc,axis=0)
            avg_inc = avg_inc
            return avg_con,avg_inc,ar_con,ar_inc

        # algorithm function returns avg_con,avg_inc,ar_con,ar_inc   
        out_final = []
        for i in range(len(channel_data.T)):
            out_final.append(algorithm(chanNames[i],scanID,trigger_channel,channel_data[:,i],period,stimTrig,clip,dispIMG))
        out_final = np.asarray(out_final,dtype="object").T
        out_final = out_final.transpose()
        return out_final

def plot_ERPs(destination_dir,data_1,data_2,latency,header,x_label,y_label,label_1,label_2,color_1,color_2,amp_range,img_name):
    """
    This plot function possesses an array of abilities
    input:    x-axis: latency
              y-axis: two 1_D arrays of data
    Functions is used to:
    1. Plot two types of erp stimulus
    2. Plot two different types of erps
    3. Plot erps for two different groups e.g., N100 for experimental and control group
    """
    fig, ax = plt.subplots()
    time = latency
    ax.plot(time, data_1,color_1, label = label_1)
    ax.plot(time, data_2,color_2, label = label_2)
    ax.set_title(header)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    fig.tight_layout()
    ax.xaxis.set_major_locator(MultipleLocator(100)) # add major ticks on x axis         
    plt.vlines(x=[0.0], ymin=amp_range,ymax=-amp_range, colors='green', ls='--',lw=2)
    ax.invert_yaxis()
    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')  
    plt.savefig(destination_dir+'/'+header+'_'+img_name+'.png',bbox_inches='tight')
    plt.show()

def spectogramPlot(data,fs,nfft,y_max,nOverlap,figsize,subTitles,title):
    #   Inputs  :   data    - 2D numpy array (d0 = samples, d1 = channels) of filtered EEG data
    #               fs      - sampling rate of hardware (defaults to config)
    #               nfft    - number of points to use in each block (defaults to config)
    #               nOverlap- number of points to overlap between blocks (defaults to config)
    #               figsize - size of figure (defaults to config)
    #               titles  - titles for each channel (defaults to config)
    y = data
    if len(y.T) % 2 != 0:
        nrows,ncols=1,int(len(y.T))
    elif len(y.T) % 2 == 0:
        nrows,ncols=2,int(len(y.T)/2)
    fig, axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=(figsize[0],figsize[1]))
    fig.suptitle(title)
    label= ["Power/Frequency"]
    for i, axs in enumerate(axs.flatten()):
        d, f, t, im = axs.specgram(data[:,i],NFFT=nfft,Fs=fs,noverlap=nOverlap)
        axs.set_title(subTitles[i])
        axs.set_ylim(0,y_max)
        #axs.set_yticks(np.arange(0,80,20))
        axs.set(xlabel='Time (s)', ylabel='Frequency (Hz)')
        axs.label_outer()
        axs
    fig.colorbar(im, ax=axs, shrink=0.9, aspect=10)

def rolling_window(data_array,timing_array,window_size,step_size):
    """
    Inputs:
    1. data_array - 1D numpy array (d0 = channels) of data
    2. timing_array - 1D numpy array (d0 = samples) of timing data
    3. len(data_array) == len(timing_array)
    4. window_size - number of samples to use in each window in seconds e.g. 1 is 1 second
    5. step_size - the step size in seconds e.g.0.5 is 0.5 seconds 

    Outputs:    
    1. data_windows - 2D numpy array (d0 = windows, d1 = window size) of data

    """
    idx_winSize = np.where(timing_array == window_size)[0][0]
    idx_stepSize = np.where(timing_array == step_size)[0][0]
    shape = (data_array.shape[0] - idx_winSize + 1, idx_winSize)
    strides = (data_array.strides[0],) + data_array.strides
    rolled = np.lib.stride_tricks.as_strided(data_array, shape=shape, strides=strides)
    return rolled[np.arange(0,shape[0],idx_stepSize)]

def slidingWindow(X, window_length, stride1):
    shape = (X.shape[0] - window_length + 1, window_length)
    strides = (X.strides[0],) + X.strides
    rolled = np.lib.stride_tricks.as_strided(X, shape=shape, strides=strides)
    return rolled[np.arange(0,shape[0],stride1)]

# Brain Vital Signs or Elemental Brain Scores functions
def normData_analysis(data,rem_outliers):

    if rem_outliers==0:
        cleaned_data = data

        def n100_amp(cleaned_data):
            n100_amp_max = np.amax(cleaned_data)
            n100_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n100_amp_best = n100_amp_max
            n100_amp_mean = np.mean(cleaned_data)
            n100_amp_std = np.std(cleaned_data)
            return n100_amp_mean,n100_amp_best,n100_amp_max,n100_amp_min,n100_amp_std

        def n100_lat(cleaned_data):
            n100_lat_max = np.amax(cleaned_data)
            n100_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n100_lat_best = n100_lat_min
            n100_lat_mean = np.mean(cleaned_data)
            n100_lat_std = np.std(cleaned_data)
            return n100_lat_mean,n100_lat_best,n100_lat_max,n100_lat_min,n100_lat_std

        def p300_amp(cleaned_data):
            p300_amp_max = np.amax(cleaned_data)
            p300_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            p300_amp_best = p300_amp_max
            p300_amp_mean = np.mean(cleaned_data)
            p300_amp_std = np.std(cleaned_data)
            return p300_amp_mean,p300_amp_best,p300_amp_max,p300_amp_min,p300_amp_std

        def p300_lat(cleaned_data):
            p300_lat_max = np.amax(cleaned_data)
            p300_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            p300_lat_best = p300_lat_min
            p300_lat_mean = np.mean(cleaned_data)
            p300_lat_std = np.std(cleaned_data)
            return p300_lat_mean,p300_lat_best,p300_lat_max,p300_lat_min,p300_lat_std
        
        def n400_amp(cleaned_data):
            n400_amp_max = np.amax(cleaned_data)
            n400_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n400_amp_best = n400_amp_max
            n400_amp_mean = np.mean(cleaned_data)
            n400_amp_std = np.std(cleaned_data)
            return n400_amp_mean,n400_amp_best,n400_amp_max,n400_amp_min,n400_amp_std

        def n400_lat(cleaned_data):
            n400_lat_max = np.amax(cleaned_data)
            n400_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n400_lat_best = n400_lat_min
            n400_lat_mean = np.mean(cleaned_data)
            n400_lat_std = np.std(cleaned_data)
            return n400_lat_mean,n400_lat_best,n400_lat_max,n400_lat_min,n400_lat_std 

    if rem_outliers==1:
        #remove outliers
        data = data.sort_values(axis=0, ascending=True)
        q1 = data.quantile(0.25)
        q3 = data.quantile(0.75)
        iqr = q3 - q1
        low = q1 - 1.5 * iqr
        high = q3 + 1.5 * iqr
        cleaned_data = data.loc[(data > low) & (data < high)]
        cleaned_data = cleaned_data.to_numpy()

        def n100_amp(cleaned_data):
            n100_amp_max = np.amax(cleaned_data)
            n100_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n100_amp_best = n100_amp_max
            n100_amp_mean = np.mean(cleaned_data)
            n100_amp_std = np.std(cleaned_data)
            return n100_amp_mean,n100_amp_best,n100_amp_max,n100_amp_min,n100_amp_std

        def n100_lat(cleaned_data):
            n100_lat_max = np.amax(cleaned_data)
            n100_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n100_lat_best = n100_lat_min
            n100_lat_mean = np.mean(cleaned_data)
            n100_lat_std = np.std(cleaned_data)
            return n100_lat_mean,n100_lat_best,n100_lat_max,n100_lat_min,n100_lat_std

        def p300_amp(cleaned_data):
            p300_amp_max = np.amax(cleaned_data)
            p300_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            p300_amp_best = p300_amp_max
            p300_amp_mean = np.mean(cleaned_data)
            p300_amp_std = np.std(cleaned_data)
            return p300_amp_mean,p300_amp_best,p300_amp_max,p300_amp_min,p300_amp_std

        def p300_lat(cleaned_data):
            p300_lat_max = np.amax(cleaned_data)
            p300_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            p300_lat_best = p300_lat_min
            p300_lat_mean = np.mean(cleaned_data)
            p300_lat_std = np.std(cleaned_data)
            return p300_lat_mean,p300_lat_best,p300_lat_max,p300_lat_min,p300_lat_std
        
        def n400_amp(cleaned_data):
            n400_amp_max = np.amax(cleaned_data)
            n400_amp_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n400_amp_best = n400_amp_max
            n400_amp_mean = np.mean(cleaned_data)
            n400_amp_std = np.std(cleaned_data)
            return n400_amp_mean,n400_amp_best,n400_amp_max,n400_amp_min,n400_amp_std

        def n400_lat(cleaned_data):
            n400_lat_max = np.amax(cleaned_data)
            n400_lat_min = np.amin(np.array(list(filter(lambda a: a != 0, cleaned_data))))
            n400_lat_best = n400_lat_min
            n400_lat_mean = np.mean(cleaned_data)
            n400_lat_std = np.std(cleaned_data)
            return n400_lat_mean,n400_lat_best,n400_lat_max,n400_lat_min,n400_lat_std 

    return cleaned_data,n100_amp,n100_lat,n400_amp,n400_lat,p300_amp,p300_lat

def ebs(parameter,mean,max,min):
    if parameter=='amplitude':
        if mean > max:
            score = 100
        if mean < min:
            score = 0
        if (mean>min and mean<max):
            best = max
            score = (round((1 - abs((mean-best)/(max-min))).item(),2))*100
    if parameter=='latency':
        if mean > max:
            score = 0
        if mean < min:
            score = 100
        if (mean>min and mean<max):
            best = min
            score = (round((1 - abs((mean-best)/(max-min))).item(),2))*100
    return score

def plots(x,y,titles,figsize,pltclr):
    x_lim = [x[0],x[-1]]
    if len(y.T) % 2 != 0:
        nrows,ncols=int(len(y.T)),1
    elif len(y.T) % 2 == 0:
        nrows,ncols=2,int(len(y.T)/2)
    fig, axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=(figsize[0],figsize[1]))
    for i, axs in enumerate(axs.flatten()):
        axs.plot(x, y[:,i], color=pltclr[i])
        axs.set_title(titles[i])
        #axs.set_ylim([np.max(y[:,i])+1000,np.min(y[:,i])-1000])
        axs.set_xlim([x_lim[0],x_lim[1]])
        axs.set(xlabel='Time (s)', ylabel='Amplitude (uV)')
        axs.tick_params(axis='both', which='major', labelsize=8)
        axs.label_outer()

def psdPlots(data,fs):
# Define window length (4 seconds)
    win = 4 * fs
    freqs,psd = signal.welch(data,fs,nperseg=win)

    # Plot the power spectrum
    sns.set(font_scale=1.2, style='white')
    plt.figure(figsize=(8, 4))
    plt.plot(freqs, psd, color='k', lw=2)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power spectral density (V^2 / Hz)')
    plt.ylim([0,0.3])
    plt.xlim([0,200])
    plt.xticks(np.arange(0,200,10))
    plt.title("Welch's periodogram")
    #plt.xlim([0, freqs.max()])
    sns.despine()

def averageERPs(device_version,chanNames,scan_IDs,dispIMG_Channel,local_path,fs,line,lowcut,highcut,stimTrig,clip,lowPassERP,label,img_name,destination_dir):
    """
    Functions generates averaged N100,P300 & N400 erps (std,dev,con & inc) from the combination of multiple eeg scans 
    Input:      1. device_version = neurocatch version
                2. scan_IDs = multiple eeg scan file names
                3. dispIMG_Channel: display channel for the ERP plot using the channel name (string) or "All" to display all channels
                4. local_path: local path to the eeg scans
    Output: averaged erps
    """
    print(label)
    def averageProcessedEEG(device_version,chanNames,scan_ID,local_path,fs,line,lowcut,highcut,stimTrig,clip,lowPassERP,dispIMG=False):
        device = importFile.neurocatch()
        fileObjects = device.init(device_version,scan_ID,local_path,dispIMG=False)
        rawEEG = fileObjects[0]
        rawEOG = fileObjects[1]
        time = fileObjects[3]
        trigOutput = fileObjects[4]
        filtering = filters()
        adaptiveFilterOutput = filtering.adaptive(rawEEG,rawEOG)
        notchFilterOutput = filtering.notch(adaptiveFilterOutput,line,fs)
        bandPassFilterOutput = filtering.butterBandPass(notchFilterOutput,lowcut,highcut,fs)
        erps = erpExtraction()
        N1P3 = erps.N100P300(chanNames,scan_ID,trigOutput,bandPassFilterOutput,time,stimTrig,clip,dispIMG=dispIMG)
        N4 = erps.N400(chanNames,scan_ID,trigOutput,bandPassFilterOutput,time,stimTrig,clip,dispIMG=dispIMG)
        N1P3_Fz,N1P3_Cz,N1P3_Pz,N4_Fz,N4_Cz,N4_Pz = N1P3[0],N1P3[1],N1P3[2],N4[0],N4[1],N4[2]
        erp_latency = np.array(np.linspace(start=-100, stop=900, num=len(N1P3_Fz[0]),dtype=object),dtype=object)
        cutoff = lowPassERP[1]
        if lowPassERP==[True,cutoff]:
            std_Fz = filtering.butter_lowpass(N1P3_Fz[0],cutoff=cutoff,fs=fs,order=4)
            dev_Fz = filtering.butter_lowpass(N1P3_Fz[1],cutoff=cutoff,fs=fs,order=4)    
            dev_Cz = filtering.butter_lowpass(N1P3_Cz[1],cutoff=cutoff,fs=fs,order=4)
            std_Cz = filtering.butter_lowpass(N1P3_Cz[0],cutoff=cutoff,fs=fs,order=4)
            dev_Pz = filtering.butter_lowpass(N1P3_Pz[1],cutoff=cutoff,fs=fs,order=4)
            std_Pz = filtering.butter_lowpass(N1P3_Pz[0],cutoff=cutoff,fs=fs,order=4)
            inc_Fz = filtering.butter_lowpass(N4_Fz[1],cutoff=cutoff,fs=fs,order=4)
            con_Fz = filtering.butter_lowpass(N4_Fz[0],cutoff=cutoff,fs=fs,order=4)
            inc_Cz = filtering.butter_lowpass(N4_Cz[1],cutoff=cutoff,fs=fs,order=4)
            con_Cz = filtering.butter_lowpass(N4_Cz[0],cutoff=cutoff,fs=fs,order=4)
            inc_Pz = filtering.butter_lowpass(N4_Pz[1],cutoff=cutoff,fs=fs,order=4)
            con_Pz = filtering.butter_lowpass(N4_Pz[0],cutoff=cutoff,fs=fs,order=4)
        elif lowPassERP==[False,False]:
            std_Fz,dev_Fz,std_Cz,dev_Cz,std_Pz,dev_Pz,inc_Fz,con_Fz,inc_Cz,con_Cz,inc_Pz,con_Pz = N1P3_Fz[0],N1P3_Fz[1],N1P3_Cz[0],N1P3_Cz[1],N1P3_Pz[0],N1P3_Pz[1],N4_Fz[0],N4_Fz[1],N4_Cz[0],N4_Cz[1],N4_Pz[0],N4_Pz[1]

        output = [std_Fz,dev_Fz,std_Cz,dev_Cz,std_Pz,dev_Pz,con_Fz,inc_Fz,con_Cz,inc_Cz,con_Pz,inc_Pz,erp_latency]
        return output
   
    device_version,chanNames,scan_IDs,local_path,fs,line,lowcut,highcut,stimTrig,clip,lowPassERP = device_version,chanNames,scan_IDs,local_path,fs,line,lowcut,highcut,stimTrig,clip,lowPassERP
    chans = dict(Fz=0,Cz=1,Pz=2,All=3)
    chan_idx = chans[dispIMG_Channel]
    if chan_idx==0:
        avgERP = []
        for i in range(len(scan_IDs)):
            avgERP.append(averageProcessedEEG(device_version=device_version,chanNames=[chanNames[0]]*len(scan_IDs),scan_ID=scan_IDs[i],local_path=local_path,fs=fs,line=line,lowcut=lowcut,highcut=highcut,stimTrig=stimTrig,clip=clip,lowPassERP=lowPassERP,dispIMG=False))
        avg_ERP = np.array(avgERP)
        avgERP = np.mean(avg_ERP,axis=0)
        plot_ERPs(destination_dir,avgERP[0],avgERP[1],avgERP[12],'N1P3_Fz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[6],avgERP[7],avgERP[12],'N4_Fz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
    if chan_idx==1:
        avgERP = []
        for i in range(len(scan_IDs)):
            avgERP.append(averageProcessedEEG(device_version=device_version,chanNames=[chanNames[1]]*len(scan_IDs),scan_ID=scan_IDs[i],local_path=local_path,fs=fs,line=line,lowcut=lowcut,highcut=highcut,stimTrig=stimTrig,clip=clip,lowPassERP=lowPassERP,dispIMG=False))
        avg_ERP = np.array(avgERP)
        avgERP = np.mean(avg_ERP,axis=0)
        plot_ERPs(destination_dir,avgERP[2],avgERP[3],avgERP[12],'N1P3_Cz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[8],avgERP[9],avgERP[12],'N4_Cz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
    if chan_idx==2:
        avgERP = []
        for i in range(len(scan_IDs)):
            avgERP.append(averageProcessedEEG(device_version=device_version,chanNames=[chanNames[2]]*len(scan_IDs),scan_ID=scan_IDs[i],local_path=local_path,fs=fs,line=line,lowcut=lowcut,highcut=highcut,stimTrig=stimTrig,clip=clip,lowPassERP=lowPassERP,dispIMG=False))
        avg_ERP = np.array(avgERP)
        avgERP = np.mean(avg_ERP,axis=0)
        plot_ERPs(destination_dir,avgERP[4],avgERP[5],avgERP[12],'N1P3_Pz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[10],avgERP[11],avgERP[12],'N4_Pz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
    if chan_idx==3:
        avgERP = []
        for i in range(len(scan_IDs)):
            avgERP.append(averageProcessedEEG(device_version=device_version,chanNames=[chanNames]*len(scan_IDs),scan_ID=scan_IDs[i],local_path=local_path,fs=fs,line=line,lowcut=lowcut,highcut=highcut,stimTrig=stimTrig,clip=clip,lowPassERP=lowPassERP,dispIMG=False))
        avg_ERP = np.array(avgERP)
        avgERP = np.mean(avg_ERP,axis=0)
        plot_ERPs(destination_dir,avgERP[0],avgERP[1],avgERP[12],'N1P3_Fz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[6],avgERP[7],avgERP[12],'N4_Fz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[2],avgERP[3],avgERP[12],'N1P3_Cz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[8],avgERP[9],avgERP[12],'N4_Cz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[4],avgERP[5],avgERP[12],'N1P3_Pz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10,img_name)
        plot_ERPs(destination_dir,avgERP[10],avgERP[11],avgERP[12],'N4_Pz','Latency (ms)','Amplitude (uV)','con','inc','b','r',10,img_name)
    return avgERP,avg_ERP

def ica(data,fs):
    """
    input: samples x channels
    output: samples x channels
    """

    #   Implement high pass filter @ 1Hz
    def icaHighpass(data,cutoff,fs):
        def params_fnc(data,cutoff,fs,order=4):
            nyq = 0.5 * fs
            normal_cutoff = cutoff / nyq
            b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
            y = signal.filtfilt(b, a, data)
            return y
        filterEEG = []
        for i in range(len(data.T)):
            filterEEG.append(params_fnc(data.T[i],cutoff,fs))
        filterEEG = np.array(filterEEG).T
        return filterEEG

    def confidenceInterval(samples):
    #   At 95% significance level, tN -1 = 2.201
        means = np.mean(samples)
        std_dev = np.std(samples)
        standard_error = std_dev/np.sqrt(len(samples))
        lower_95_perc_bound = means - 2.201*standard_error
        upper_95_perc_bound = means + 2.201*standard_error
        return upper_95_perc_bound

    def setZeros(data,index):
        def params(data):
            return np.zeros(len(data))
        zeros = []
        for i in range(len(index)):
            zeros.append(params(data.T[index[i]]))
        zeros = np.array(zeros)
        return zeros

    hpEEG = icaHighpass(data,cutoff=1,fs=fs) 

    #   Computing ICA components
    ica = FastICA(n_components=3, random_state=0, tol=0.0001)
    comps = ica.fit_transform(hpEEG)
    comps_1 = comps[:,0]
    comps_2 = comps[:,1]
    comps_3 = comps[:,2]

    #   Computing kurtosis of ICA weights
    comps_1_kurtosis = kurtosis(comps_1)
    comps_2_kurtosis = kurtosis(comps_2)
    comps_3_kurtosis = kurtosis(comps_3)
    comps_kurtosis = np.array([comps_1_kurtosis,comps_2_kurtosis,comps_3_kurtosis])

    #   Computing skewness of ICA weights
    comps_1_skew = skew(comps_1)
    comps_2_skew = skew(comps_2)
    comps_3_skew = skew(comps_3)
    comps_skew = np.array([comps_1_skew,comps_2_skew,comps_3_skew])

    #   Computing sample entropy of ICA weights
    import antropy as ant
    comps_1_sampEN = ant.sample_entropy(comps_1)
    comps_2_sampEN = ant.sample_entropy(comps_2)
    comps_3_sampEN = ant.sample_entropy(comps_3)
    comps_sampEN = np.array([comps_1_sampEN,comps_2_sampEN,comps_3_sampEN])

    #   Computing CI on to set threshold
    threshold_kurt = confidenceInterval(comps_kurtosis)
    threshold_skew = confidenceInterval(comps_skew)
    threshold_sampEN = confidenceInterval(comps_sampEN)

    "compare threshold with extracted parameter values"
    #   Extract epochs
    bool_ArtfCompsKurt = [comps_kurtosis>threshold_kurt]
    idx_ArtfCompsKurt = np.asarray(np.where(bool_ArtfCompsKurt[0]==True))
    bool_ArtfCompsSkew = [comps_skew>threshold_skew]
    idx_ArtfCompsSkew = np.asarray(np.where(bool_ArtfCompsSkew[0]==True))
    bool_ArtfCompsSampEN = [comps_sampEN>threshold_sampEN]
    idx_ArtfCompsSampEN = np.asarray(np.where(bool_ArtfCompsSampEN[0]==True))

    #   Merge index of components detected as artifacts by kurtosis, skewness, and sample entropy
    idx_artf_comps = np.concatenate((idx_ArtfCompsKurt,idx_ArtfCompsSkew,idx_ArtfCompsSampEN),axis=1)
    idx_artf_comps = np.unique(idx_artf_comps)

    "Component identified as artifact is converted to arrays of zeros"
    rejected_comps = setZeros(comps,idx_artf_comps)


    "Return zero-ed ICs into the original windows per ICs"
    for i in range(len(idx_artf_comps)):
        idx_rejected_comps = np.arange(len(rejected_comps))
        comps.T[idx_artf_comps[i]] = rejected_comps[idx_rejected_comps[i]]


    "Recover clean signal from clean ICs"
    restored = ica.inverse_transform(comps)
    return restored