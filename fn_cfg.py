# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:00:13 2021

@author: oseho
"""


from df_lib import*
from params import*


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
        def init(self,version,filename,localPath):
            if version == 1.0:
                data = filename
                localPath = localPath.replace(os.sep, '/')   
                localPath = localPath + '/'  
                path = localPath+data
                os.chdir(path)
                for file in (os.listdir(path))[0:2]:
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
                for i in metadata:
                    print(i)
                print(metadata['version'])
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
            print(filenames)
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
                print('channel P3 utilized')
            else:
                eogChan1 = np.zeros(len(fz.T))
                eogChan1 = eogChan1.reshape(1,len(eogChan1))

            if p4 < 501:
                eogChans = gtec['eogChans']
                eog20 = eogChans[1]
                eogChan2 = data1[:,eog20]
                eogChan2 = eogChan2.reshape(1,dataRows)
                print('channel P4 utilized')
            else:
                eogChan2 = np.zeros(len(fz.T))
                eogChan2 = eogChan2.reshape(1,len(eogChan2))

            if p07 < 501:
                eogChans = gtec['eogChans']
                eog30 = eogChans[2]
                eogChan3 = data1[:,eog30]
                eogChan3 = eogChan3.reshape(1,dataRows)
                print('channel P07 utilized')
            else:
                eogChan3 = np.zeros(len(fz.T))
                eogChan3 = eogChan3.reshape(1,len(eogChan3))

            if p08 < 501:
                eogChans = gtec['eogChans']
                eog40 = eogChans[3]
                eogChan4 = data1[:,eog40]
                eogChan4 = eogChan4.reshape(1,dataRows)
                print('channel P08 utilized')
            else:
                eogChan4 = np.zeros(len(fz.T))
                eogChan4 = eogChan4.reshape(1,len(eogChan4))

            if oz < 501:
                eogChans = gtec['eogChans']
                eog50 = eogChans[4]
                eogChan5 = data1[:,eog50]
                eogChan5 = eogChan5.reshape(1,dataRows)
                print('channel 0Z utilized')
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
                rawEEGEOG = np.concatenate((rawEEG),axis=1)
                print('data contains Fz, Cz, Pz & no EOG channels')
            
            elif len(rawData.T)==5:
                csm_ntrig = 4
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                print('data contains Fz, Cz, Pz & one EOG channel')

            elif len(rawData.T)==6:
                csm_ntrig = 5
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                print('data contains Fz, Cz, Pz & two EOG channels')

            elif len(rawData.T)==7:
                csm_ntrig = 6
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                print('data contains Fz, Cz, Pz & three EOG channels')

            elif len(rawData.T)==8:
                csm_ntrig = 7
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3,csm_eog4]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                print('data contains Fz, Cz, Pz & four EOG channels')

            elif len(rawData.T)==9:
                csm_ntrig = 8
                rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_eog5,csm_ntrig]]  
                rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
                rawEOG = rawData[:,[csm_eog1,csm_eog2,csm_eog3,csm_eog4,csm_eog5]]
                rawEEGEOG = np.concatenate((rawEEG,rawEOG),axis=1)
                print('data contains Fz, Cz, Pz & five EOG channels')

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
    def notch(self,data,lines,fs,Q=30):
        """
           Inputs  :   data    - 2D numpy array (d0 = samples, d1 = channels) of unfiltered EEG data
                       cut     - frequency to be notched (defaults to config)
                       fs      - sampling rate of hardware (defaults to config)
                       Q       - Quality Factor (defaults to 30) that characterizes notch filter -3 dB bandwidth bw relative to its center frequency, Q = w0/bw.   
           Output  :   y     - 2D numpy array (d0 = samples, d1 = channels) of notch-filtered EEG data
           NOTES   :   
           Todo    : report testing filter characteristics
        """
        def initialNotch(data,lines,fs,Q=30):
            cut = line
            w0 = cut/(fs/2)
            b, a = signal.iirnotch(w0, Q)
            y = signal.filtfilt(b, a, data, axis=0)
            return y
        y_1 = initialNotch(data,lines[0],fs,Q)
        y_2 = initialNotch(y_1,lines[1],fs,Q)
        y_3 = initialNotch(y_2,lines[2],fs,Q)
        y_4 = initialNotch(y_3,lines[3],fs,Q)
        return y_4

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
    def N100P300(self,trigger_channel,eegData,period,stimTrig,clip):
        """
          Inputs: trigger channels, bandpass filtered data for all channels, time period,
                  stimTrig, clip value
        """
        trigger_channel,channel_data,period,stimTrig,clip = trigger_channel,eegData,period,stimTrig,clip
        def algorithm(trigger_channel,channel_data,period,stimTrig,clip):
            trigger_data = rising_edge(trigger_channel)
            trigCol = trigger_data
            avg = channel_data 
            Ts = period
            # STANDARD TONE [1]: extract the time points where stimuli 1 exists
            std = stimTrig['std']
            std = std[0]
            
            no_ones = np.count_nonzero(trigCol==1)
            print("number of std tone event codes:",no_ones)
            
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
                l = endVals
                max_len = max([len(arr) for arr in l])
                padded = np.array([np.lib.pad(arr, (0, max_len - len(arr)), 'constant', constant_values=0) for arr in l])
                    # appropriate end columns
                endCols = padded
                
                # merge the different sections of the epoch to form the epochs
                trig_stdpad = np.concatenate((startCols,desiredCols,endCols),axis=1)
            
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
                l_eeg = eeg_endVals
                eeg_max_len = max([len(arr) for arr in l_eeg])
                eeg_padded = np.array([np.lib.pad(arr, (0, eeg_max_len - len(arr)), 'constant', constant_values=0) for arr in l_eeg])
                    # appropriate end columns
                eeg_endCols = eeg_padded
                
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
                print("notice! epochs lost for std tone:","{:.2%}".format((int(dif))/100))
            elif len(ar_std) < len(bc_std):
                print("callback! epochs lost for std tone:","{:.2%}".format((int(dif))/100))
            
                
            # averaging
            avg_std = np.mean(ar_std,axis=0)

            #%%
            # DEVIANT TONE [2]: extract the time points where stimuli 2 exists
            dev = stimTrig['dev']
            dev = dev[0]
            
            no_twos = np.count_nonzero(trigCol==2)
            print("number of dev tone event codes:",no_twos)
            
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
                l = endVals
                max_len = max([len(arr) for arr in l])
                padded = np.array([np.lib.pad(arr, (0, max_len - len(arr)), 'constant', constant_values=0) for arr in l])
                    # appropriate end columns
                endCols = padded
                
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
                l_eeg = eeg_endVals
                eeg_max_len = max([len(arr) for arr in l_eeg])
                eeg_padded = np.array([np.lib.pad(arr, (0, eeg_max_len - len(arr)), 'constant', constant_values=0) for arr in l_eeg])
                    # appropriate end columns
                eeg_endCols = eeg_padded
                
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
                print("notice! epochs lost for dev tone:","{:.2%}".format((int(dif))/100))
            elif len(ar_dev) < len(bc_dev):
                print("callback! epochs lost for dev tone:","{:.2%}".format((int(dif))/100))
            
            # averaging
            avg_dev = np.mean(ar_dev,axis=0)
            return avg_std,avg_dev,ar_std,ar_dev
        # algorithm function returns avg_std,avg_dev,ar_std,ar_dev
        out_final = []
        for i in range(len(channel_data.T)):
            out_final.append(algorithm(trigger_channel,channel_data[:,i],period,stimTrig,clip))
        out_final = np.asarray(out_final).T
        out_final = out_final.transpose()
        return out_final

    def N400(self,trigger_channel,eegData,period,stimTrig,clip):
        """
          Inputs: trigger channels, bandpass filtered data for all channels, time period,
                  stimTrig, clip value
        """
        trigger_channel,channel_data,period,stimTrig,clip = trigger_channel,eegData,period,stimTrig,clip
        def algorithm(trigger_channel,channel_data,period,stimTrig,clip):
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
            print("number of con word event codes:",no_cons)
            
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
                l = endVals
                max_len = max([len(arr) for arr in l])
                padded = np.array([np.lib.pad(arr, (0, max_len - len(arr)), 'constant', constant_values=0) for arr in l])
                    # appropriate end columns
                endCols = padded
                
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
                l_eeg = eeg_endVals
                eeg_max_len = max([len(arr) for arr in l_eeg])
                eeg_padded = np.array([np.lib.pad(arr, (0, eeg_max_len - len(arr)), 'constant', constant_values=0) for arr in l_eeg])
                    # appropriate end columns
                eeg_endCols = eeg_padded
                
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
                print("notice! epochs lost for con word:","{:.2%}".format((int(dif))/100))
            elif len(ar_con) < len(bc_con):
                print("callback! epochs lost for con word:","{:.2%}".format((int(dif))/100))
            
            # averaging
            avg_con = np.mean(ar_con,axis=0)
            
            
            # %%
            # incongruent word [4,7]: extract the time points where stimuli 1 exists
            inc = stimTrig['inc']
            inc = inc[0:2]
            inc = (np.array([inc])).T
            
            no_fives = np.count_nonzero(trigCol==5)
            no_eights = np.count_nonzero(trigCol==8)
            no_incs = no_fives + no_eights
            print("number of inc word event codes:",no_incs)
            
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
                l_eeg = eeg_endVals
                eeg_max_len = max([len(arr) for arr in l_eeg])
                eeg_padded = np.array([np.lib.pad(arr, (0, eeg_max_len - len(arr)), 'constant', constant_values=0) for arr in l_eeg])
                    # appropriate end columns
                eeg_endCols = eeg_padded
                
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
                print("notice! epochs lost for inc word:","{:.2%}".format((int(dif))/100))
            elif len(ar_inc) < len(bc_inc):
                print("callback! epochs lost for inc word:","{:.2%}".format((int(dif))/100))
            # averaging
            avg_inc = np.mean(ar_inc,axis=0)
            return avg_con,avg_inc,ar_con,ar_inc

        # algorithm function returns avg_con,avg_inc,ar_con,ar_inc   
        out_final = []
        for i in range(len(channel_data.T)):
            out_final.append(algorithm(trigger_channel,channel_data[:,i],period,stimTrig,clip))
        out_final = np.asarray(out_final).T
        out_final = out_final.transpose()
        return out_final

def plot_ERPs(data_1,data_2,latency,header,x_label,y_label,label_1,label_2,color_1,color_2,amp_range):
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
    #ax.yaxis.grid()
    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')  
    plt.show()

def neurocatchPipeline(filename,localPath,line,fs,Q,stimTrig):
    # combines preprocessing and post processing
    # processes erps for each channel
    raw_data = importFile('neurocatch_1_2')
    raw_data = raw_data.neurocatch_1_2(filename,localPath)
    ts = raw_data[2]
    adp_data = adaptive_filter(raw_data[1],raw_data[4],ts,plot='false')
    not_data = notch_filter('neurocatch_1_2')
    not_data = not_data.neurocatch_1_2(adp_data,line,fs,Q)
    bpData = butter_bandpass_filter('neurocatch_1_2')
    bp_data = bpData.neurocatch(not_data)
    trig_chan = raw_data[3]
    trig_xform = rising_edge(trig_chan)
    # fz channel
    fz_tones = tones(trig_xform,bp_data[:,0],ts,stimTrig)
    fz_words = words(trig_xform,bp_data[:,0],ts,stimTrig)
    # cz channel
    cz_tones = tones(trig_xform,bp_data[:,1],ts,stimTrig)
    cz_words = words(trig_xform,bp_data[:,1],ts,stimTrig)
    # pz channel
    pz_tones = tones(trig_xform,bp_data[:,2],ts,stimTrig)
    pz_words = words(trig_xform,bp_data[:,2],ts,stimTrig)
    erp_latency = np.array(np.linspace(start=-100, stop=900, num=len(fz_tones[0])))
    # tones[0] = std; tones[1] = dev; words[0] = con; words[1] = inc
    return fz_tones[0],fz_tones[1],fz_words[0],fz_words[1],cz_tones[0],cz_tones[1],cz_words[0],cz_words[1],pz_tones[0],pz_tones[1],pz_words[0],pz_words[1],erp_latency,fz_tones[2],fz_tones[3]

def fz_gndavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[0]]))
        zeros_dev = np.vstack(([zeros_dev, erp[1]]))
        zeros_con = np.vstack(([zeros_con, erp[2]]))
        zeros_inc = np.vstack(([zeros_inc, erp[3]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        avg_std = np.nanmean(scans_std,axis = 1)
        return scans_std,avg_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        avg_dev = np.nanmean(scans_dev,axis = 1)
        return scans_dev,avg_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        avg_con = np.nanmean(scans_con,axis = 1)
        return scans_con,avg_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        avg_inc = np.nanmean(scans_inc,axis = 1)
        return scans_inc,avg_inc,erp[12]
    
    return std,dev,con,inc

def cz_gndavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[4]]))
        zeros_dev = np.vstack(([zeros_dev, erp[5]]))
        zeros_con = np.vstack(([zeros_con, erp[6]]))
        zeros_inc = np.vstack(([zeros_inc, erp[7]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        avg_std = np.nanmean(scans_std,axis = 1)
        return scans_std,avg_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        avg_dev = np.nanmean(scans_dev,axis = 1)
        return scans_dev,avg_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        avg_con = np.nanmean(scans_con,axis = 1)
        return scans_con,avg_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        avg_inc = np.nanmean(scans_inc,axis = 1)
        return scans_inc,avg_inc,erp[12]
    
    return std,dev,con,inc

def pz_gndavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[8]]))
        zeros_dev = np.vstack(([zeros_dev, erp[9]]))
        zeros_con = np.vstack(([zeros_con, erp[10]]))
        zeros_inc = np.vstack(([zeros_inc, erp[11]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        avg_std = np.nanmean(scans_std,axis = 1)
        return scans_std,avg_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        avg_dev = np.nanmean(scans_dev,axis = 1)
        return scans_dev,avg_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        avg_con = np.nanmean(scans_con,axis = 1)
        return scans_con,avg_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        avg_inc = np.nanmean(scans_inc,axis = 1)
        return scans_inc,avg_inc,erp[12]
    
    return std,dev,con,inc

def fz_sngavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[0]]))
        zeros_dev = np.vstack(([zeros_dev, erp[1]]))
        zeros_con = np.vstack(([zeros_con, erp[2]]))
        zeros_inc = np.vstack(([zeros_inc, erp[3]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        return scans_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        return scans_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        return scans_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        return scans_inc,erp[12]
    
    return std,dev,con,inc

def cz_sngavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[4]]))
        zeros_dev = np.vstack(([zeros_dev, erp[5]]))
        zeros_con = np.vstack(([zeros_con, erp[6]]))
        zeros_inc = np.vstack(([zeros_inc, erp[7]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        return scans_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        return scans_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        return scans_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        return scans_inc,erp[12]
    
    return std,dev,con,inc

def pz_sngavg(scan_id,no_scans,len_epoch,gtec):
    scans = scan_id
    zeros_std = np.zeros(len_epoch)
    zeros_dev = np.zeros(len_epoch)
    zeros_con = np.zeros(len_epoch)
    zeros_inc = np.zeros(len_epoch)
    for i in range(no_scans):
        erp = pre_post_process(scans[i],gtec)
        zeros_std = np.vstack(([zeros_std, erp[8]]))
        zeros_dev = np.vstack(([zeros_dev, erp[9]]))
        zeros_con = np.vstack(([zeros_con, erp[10]]))
        zeros_inc = np.vstack(([zeros_inc, erp[11]]))

    def std():
        scans_std = np.delete(zeros_std,0,axis=0).T
        return scans_std,erp[12]

    def dev():
        scans_dev = np.delete(zeros_dev,0,axis=0).T
        return scans_dev,erp[12]

    def con():
        scans_con = np.delete(zeros_con,0,axis=0).T
        return scans_con,erp[12]

    def inc():
        scans_inc = np.delete(zeros_inc,0,axis=0).T
        return scans_inc,erp[12]
    
    return std,dev,con,inc

def ptp_erpscan(peak_val,erp_data,subjs_data):
    p2p = peaktopeak(erp_data.T)
    result_1 = (np.where(p2p > peak_val))[0]
    acpt_erp = np.delete(erp_data.T,(result_1),axis = 0)
    result_2 = [~np.isnan(acpt_erp).any(axis=1)]
    acpt_erp = acpt_erp[result_2]
    acpt_subjs = np.delete(subjs_data,(result_1),axis = 0)
    acpt_subjs = acpt_subjs[result_2]
    return acpt_erp,acpt_subjs

def spectogramPlot(data,fs,nfft,nOverlap,figsize,subTitles,title):
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
        axs.set_ylim(0,np.amax(f))
        axs.set_yticks(np.arange(0,np.amax(f),20))
        axs.set(xlabel='Time (s)', ylabel='Frequency (Hz)')
        axs.label_outer()
        axs
    fig.colorbar(im, ax=axs, shrink=0.9, aspect=10)

# signal quality evaluating functions
def rolling_window(array, window_size,freq):
    shape = (array.shape[0] - window_size + 1, window_size)
    strides = (array.strides[0],) + array.strides
    rolled = np.lib.stride_tricks.as_strided(array, shape=shape, strides=strides)
    return rolled[np.arange(0,shape[0],freq)]

def general_amplitude(data,x,y):    
    def genAmp(data):
        bins = np.arange(-1000,1000,10)
        counts,bin_edges = np.histogram(data, bins)
        #    plt.hist(counts,bins)
        #    plt.show()
        # index of -100uV
        idx_minUV = int((np.where(bins == -100))[0])
        # index of 100uV
        idx_maxUV = int((np.where(bins == 100))[0])
        # extract counts between -100uV to 100uV
        nrmCounts = counts[idx_minUV:idx_maxUV]
        # extract the highest count between counts of -100uV to 100uV
        max_counts = np.amax(nrmCounts)
        # locate index of highest count
        result = np.where(nrmCounts==max_counts)
        result = result[0]
        if result.size == 1:
            idx_maxCount = np.asscalar(result)
            # check for increment in first half of bin distribution between -100uV to 100uV
            # fH = first half of bin distribution
            fH = nrmCounts[0:(idx_maxCount+1)]
            zero = np.array([0])
            concat_fH = np.hstack((zero,fH))
            diff_fH = np.diff(concat_fH)
            # check for peak change between count of bins 
            pc_fH = diff_fH > 0
            # extract the elements with "true" peak changes
            idx_true1 = np.where(pc_fH == True)
            # use true indices to extract counts from fH
            nbins_fH = fH[idx_true1]
            # count the total number of bins in first half
            idx_mc = np.where(counts==max_counts)
            idx_mc = idx_mc[0]
            idx_mc = idx_mc[0]
            tbins_fH = counts[0:idx_mc+1]
            # evaluate score for first half of bin
            s1 = np.sum(nbins_fH)/np.sum(tbins_fH)
            # incase it run into "Nan" scores
            if np.isnan(s1):
                s1 = np.nan_to_num(s1)
            
            # check for increment in second half of bin distribution between -100uV to 100uV
            # sH = second half of bin distribution
            sH = nrmCounts[idx_maxCount:len(nrmCounts)]
            zero = np.array([0])
            concat_sH = np.hstack((sH,zero))
            diff_sH = np.diff(concat_sH)
            # check for peak change between count of bins 
            pc_sH = diff_sH < 0
            # extract the elements with "true" peak changes
            idx_true2 = np.where(pc_sH == True)
            # use true indices to extract normal bin counts from fH
            nbins_sH = sH[idx_true2]
            # count the total number of bins in first half
            idx_mc = np.where(counts==max_counts)
            idx_mc = idx_mc[0]
            idx_mc = idx_mc[0]
            tbins_sH = counts[idx_mc:len(counts)]
            # evaluate score for first half of bin
            s2 = np.sum(nbins_sH)/np.sum(tbins_sH)
            # incase it run into "Nan" scores
            if np.isnan(s2):
                s2 = np.nan_to_num(s2)
            s3 = int(((s1+s2)/2)*100)
            
        elif result.size > 1:
            # in case of more than 1 highest count e.g., three, we pick the first highest count
            idx_maxCount = result[0]
            # check for increment in first half of bin distribution between -100uV to 100uV
            # fH = first half of bin distribution
            fH = nrmCounts[0:(idx_maxCount+1)]
            zero = np.array([0])
            concat_fH = np.hstack((zero,fH))
            diff_fH = np.diff(concat_fH)
            # check for peak change between count of bins 
            pc_fH = diff_fH > 0
            # extract the elements with "true" peak changes
            idx_true1 = np.where(pc_fH == True)
            # use true indices to extract counts from fH
            nbins_fH = fH[idx_true1]
            # count the total number of bins in first half
            idx_mc = np.where(counts==max_counts)
            idx_mc = idx_mc[0]
            idx_mc = idx_mc[0]
            tbins_fH = counts[0:idx_mc+1]
            # evaluate score for first half of bin
            s1 = np.sum(nbins_fH)/np.sum(tbins_fH)
            # incase it run into "Nan" scores
            if np.isnan(s1):
                s1 = np.nan_to_num(s1)
            
            # check for increment in second half of bin distribution between -100uV to 100uV
            # we pick the third (or last) highest count to create the second half
            idx_maxCount = result[(result.size-1)]
            # sH = second half of bin distribution
            sH = nrmCounts[idx_maxCount:len(nrmCounts)]
            zero = np.array([0])
            concat_sH = np.hstack((sH,zero))
            diff_sH = np.diff(concat_sH)
            # check for peak change between count of bins 
            pc_sH = diff_sH < 0
            # extract the elements with "true" peak changes
            idx_true2 = np.where(pc_sH == True)
            # use true indices to extract normal bin counts from fH
            nbins_sH = sH[idx_true2]
            # count the total number of bins in first half
            idx_mc = np.where(counts==max_counts)
            idx_mc = idx_mc[0]
            idx_mc = idx_mc[(idx_mc.size-1)]
            tbins_sH = counts[idx_mc:len(counts)]
            # evaluate score for first half of bin
            s2 = np.sum(nbins_sH)/np.sum(tbins_sH)
            # incase it run into "Nan" scores
            if np.isnan(s2):
                s2 = np.nan_to_num(s2)
            s3 = int(((s1+s2)/2)*100)
        return s3

    window = rolling_window(data,x,y)
    
    com_arr = np.array([])
    for x in range(len(window)): 
        sub_arr = genAmp(window[x,:])
        com_arr = np.hstack([com_arr, sub_arr])
    com_arr = com_arr.reshape((len(com_arr)),1) 
    
    rank = (np.where(com_arr < 100))[0]
    win_score = np.where(com_arr == 100, 0, com_arr)
    win_score = np.where(win_score > 0, 1, win_score)
    win_score = win_score.reshape(1,len(win_score))
    clean_win = np.count_nonzero(win_score==0)
    score1 = int((clean_win/len(win_score.T))*100)
    print("score 1:",score1)
    return score1,win_score

def amplitude_spectrum(v,w,x,y):    
    data = v
    fs = w
    def ampspect(data,fs):
        # fs = 500
        win = 4 * fs
        freq, psd= signal.welch(data, fs, nperseg=win)
        low, high = 1, 4
        idx_delta = np.logical_and(freq >= low, freq <= high)
        result = np.where(idx_delta == True)
        idx = result[0]
        deltaPower = psd[idx]

        freq, psd= signal.welch(data, fs, nperseg=win)
        low, high = 4, 8
        idx_theta = np.logical_and(freq >= low, freq <= high)
        result = np.where(idx_theta == True)
        idx = result[0]
        thetaPower = psd[idx]

        freq, psd= signal.welch(data, fs, nperseg=win)
        low, high = 8, 12
        idx_alpha = np.logical_and(freq >= low, freq <= high)
        result = np.where(idx_alpha == True)
        idx = result[0]
        alphaPower = psd[idx]

        freq, psd= signal.welch(data, fs, nperseg=win)
        low, high = 12, 30
        idx_beta = np.logical_and(freq >= low, freq <= high)
        result = np.where(idx_beta == True)
        idx = result[0]
        betaPower = psd[idx]

        freq, psd= signal.welch(data, fs, nperseg=win)
        low, high = 30, 50
        idx_gamma = np.logical_and(freq >= low, freq <= high)
        result = np.where(idx_gamma == True)
        idx = result[0]
        gammaPower = psd[idx]
        
        bw_bands = np.concatenate((deltaPower, thetaPower, alphaPower, betaPower, gammaPower),axis=0)
        avg_spectrum = np.mean(bw_bands)
        return avg_spectrum
    
    window = rolling_window(data,x,y)
    
    com_arr = np.array([])
    for x in range(len(window)): 
        sub_arr = ampspect(window[x,:],fs)
        com_arr = np.hstack([com_arr, sub_arr])
    com_arr = com_arr.reshape((len(com_arr)),1) 

    amp_spec = com_arr
    mean = np.mean(data)
    std = np.std(data)
    x = amp_spec

    win_score4 = np.array([])

    for i in range(len(x)):
        if ((mean) <= (x[i,:]) <= (1*std)) or ((mean) >= (x[i,:]) >= (-1*std)):
            k4 = 0
        elif ((x[i,:]) >= (1*std)) or ((x[i,:]) <= (-1*std)):
            k4 = 1
        win_score4 = np.hstack([win_score4,k4])

    count_zeros = np.count_nonzero(win_score4==0)
    s2 = int(((count_zeros)/(len(win_score4)))*100)
    print("score 2:",s2)
    win_score4 = win_score4.reshape(1,len(win_score4))
    return s2,win_score4

def maximum_gradient(data,x,y):     
    def MG(data):
        res = [data[i + 1] - data[i] for i in range(len(data)-1)]
        res = np.asarray(res)
        maxGrad = np.max(res)
        return maxGrad
    
    window = rolling_window(data,x,y)
    
    com_arr = np.array([])
    for x in range(len(window)): 
        sub_arr = MG(window[x,:])
        com_arr = np.hstack([com_arr, sub_arr])
    com_arr = com_arr.reshape((len(com_arr)),1) 

    maxGrad = com_arr
    mean = np.mean(data)
    std = np.std(data)
    x = maxGrad

    win_score6 = np.array([])

    for i in range(len(x)):
        if ((mean) <= (x[i,:]) <= (1*std)) or ((mean) >= (x[i,:]) >= (-1*std)):
            k6 = 1
        elif ((x[i,:]) >= (1*std)) or ((x[i,:]) <= (-1*std)):
            k6 = 0
        win_score6 = np.hstack([win_score6,k6])

    count_zeros = np.count_nonzero(win_score6==0)
    s3 = int(((count_zeros)/(len(win_score6)))*100)
    print("score 3:",s3)
    win_score6 = win_score6.reshape(1,len(win_score6))
    return s3,win_score6

def total_quality(x,y,z):
    stackedMetrics = (np.concatenate((x,y,z),axis = 0))
    sumStackedMetrics = np.mean(stackedMetrics,axis = 0)
    sumStackedMetrics = sumStackedMetrics.reshape((len(sumStackedMetrics)),1)
    sumStackedMetrics = sumStackedMetrics.T   # metrics variable hold the 0's, 1's or 2's
    total_score = np.count_nonzero(sumStackedMetrics==0)
    eegQuality = int(((total_score)/(len(stackedMetrics.T)))*100)
    print("amount of clean signal:","{:.2%}".format((int(eegQuality))/100))
    return eegQuality

def print_array(arr):
    """
    prints a 2-D numpy array in a nicer format
    """
    for a in arr:
        for elem in a:
            print("{}".format(elem).rjust(3), end="    ")
        print(end="\n")

def sq_chan(data):
    idx_folders = data
    info1 = idx_folders[0:7]
    info2 = idx_folders[7:20]
        # %% extract correct eog position
    def listToString(s):  
           # initialize an empty string 
           str1 = ""  
           # traverse in the string   
           for ele in s:  
               str1 += ele   
       # return string   
           return str1
       # Driver code     
    
    # EOG channel extraction using subject metadata i.e., eogChan_locator algorithm
    t = r"/Users/oseho/sfuvault/Documents/brainNet/projects/laurel place/dev/dataset/"
    p = info1 # info1
    u = info2 # info 2
    txt = [t,p,u] 
    txt = (listToString(txt))
    txt = [txt]
    etxt = txt[0]
    
    path = etxt
    
    # Change the directory
    os.chdir(path)
    
    # iterate through all file
    for file in os.listdir():
        # Check whether file is in text format or not
        if file.endswith(".txt"):
            file_path = f"{path}\{file}"
            
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
       
    p3 = lines[9] # extract line
    if len(p3)<=15:
        p3 = two_digits(p3)
    else:
        p3 = three_digits(p3)
    
    p4 = lines[11]
    if len(p4)<=15:
        p4 = two_digits(p4)
    else:
        p4 = three_digits(p4)
    
    
    p07 = lines[12]
    if len(p07)<=15:
        p07 = two_digits(p07)
    else:
        p07 = three_digits(p07)
    
    
    p08 = lines[13]
    if len(p08)<=15:
        p08 = two_digits(p08)
    else:
        p08 = three_digits(p08)
    
    oz = lines[14]
    if len(oz)<=15:
        oz = two_digits(oz)
    else:
        oz = three_digits(oz)
    
    
    # %% import raw file 1
    t = "/Users/oseho/sfuvault/Documents/brainNet/projects/laurel place/dev/dataset/"
    v = '.bin'
    y = '/'
    s = [t,p,u,y,p,u,v] 
    r = (listToString(s))  
    r = [r]
    filenames = glob.glob(r[0])
    print(filenames)
    data = [np.fromfile(f, dtype=np.float32) for f in filenames]
    data1 = data[0]
    
    dataCols = gtec['dataCols']
    dataRows = int(len(data1)/dataCols)                                # 14 columns or channels
    data1 = data1.reshape(dataRows, dataCols)
    
    # %% configure eeg channels
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
    else:
        eogChan1 = np.zeros(len(fz.T))
        eogChan1 = eogChan1.reshape(1,len(eogChan1))
    
    if p4 < 501:
        eogChans = gtec['eogChans']
        eog20 = eogChans[1]
        eogChan2 = data1[:,eog20]
        eogChan2 = eogChan2.reshape(1,dataRows)
    else:
        eogChan2 = np.zeros(len(fz.T))
        eogChan2 = eogChan2.reshape(1,len(eogChan2))
    
    if p07 < 501:
        eogChans = gtec['eogChans']
        eog30 = eogChans[2]
        eogChan3 = data1[:,eog30]
        eogChan3 = eogChan3.reshape(1,dataRows)
    else:
        eogChan3 = np.zeros(len(fz.T))
        eogChan3 = eogChan3.reshape(1,len(eogChan3))
    
    if p08 < 501:
        eogChans = gtec['eogChans']
        eog40 = eogChans[3]
        eogChan4 = data1[:,eog40]
        eogChan4 = eogChan4.reshape(1,dataRows)
    else:
        eogChan4 = np.zeros(len(fz.T))
        eogChan4 = eogChan4.reshape(1,len(eogChan4))
    
    if oz < 501:
        eogChans = gtec['eogChans']
        eog50 = eogChans[4]
        eogChan5 = data1[:,eog50]
        eogChan5 = eogChan5.reshape(1,dataRows)
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
    csm = dict(fz=[0], cz=[1], pz=[2], eog1=[3], eog2=[4], ntrig=[5])
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
    csm_ntrig = csm['ntrig']
    csm_ntrig = csm_ntrig[0]


    # condition 1 assesses signal quality of raw eeg data without eog channels
    if len(rawData.T)==4:
        csm_ntrig = 3
        rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_ntrig]]  
        rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              
        print("EOG data not available for this subject")
        eegQuality1 = 0
        eegQuality2 = 0
        eegQuality3 = 0
        eegQuality = (np.array((eegQuality1,eegQuality2,eegQuality3)))


    # condition 2 assesses signal quality of raw eeg data with just one eog channel
    elif len(rawData.T)==5:
        csm_ntrig = 4
        rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_ntrig]]  
        rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              # rawData4 is rawData3 without trigger channel (index 5)
        rawEOG = rawData[:,[csm_eog1]]

        # time period of scan
        fs = gtec['fs']
        sfreq = fs
        dt = 1/sfreq
        stop = dataRows/sfreq 
        Ts = np.arange(0,stop,dt) 
        Ts = Ts.reshape((len(Ts)),1)
    
        # %% Adaptive filter
        eegData = rawEEG
        eogData = rawEOG
        forgetF=0.995
        nKernel=5
        startSample=0
        p = False
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
        
        # %% EEG Signal Quality Metrics: first time domain metric using score 1 is based on the rms  
        fz_Q1 = general_amplitude(cleanData[:,0],2000,500)
        cz_Q1 = general_amplitude(cleanData[:,1],2000,500)
        pz_Q1 = general_amplitude(cleanData[:,2],2000,500) # (len(x.T)) is equal to the number of windows
        
        # %%  EEG Signal Quality Metrics: second time domain metric using score 3  : zero crossing rate
        fz_Q2 = amplitude_spectrum(cleanData[:,0],fs,2000,500)
        cz_Q2 = amplitude_spectrum(cleanData[:,1],fs,2000,500)
        pz_Q2 = amplitude_spectrum(cleanData[:,2],fs,2000,500)
        
        # %%  EEG Signal Quality Metrics: score 6: Kurtosis
        fz_Q3 = maximum_gradient(cleanData[:,0],2000,500)
        cz_Q3 = maximum_gradient(cleanData[:,1],2000,500)
        pz_Q3 = maximum_gradient(cleanData[:,2],2000,500)

        # Total score for fz channel
        fz_Quality = total_quality(fz_Q1[1], fz_Q2[1], fz_Q3[1])
        cz_Quality = total_quality(cz_Q1[1], cz_Q2[1], cz_Q3[1])
        pz_Quality = total_quality(pz_Q1[1], pz_Q2[1], pz_Q3[1])
        eegQuality = (np.array((fz_Quality,cz_Quality,pz_Quality)))


    # condition 3 assesses signal quality of raw eeg data with just two eog channel
    elif len(rawData.T)==6:
        rawData = rawData[:,[csm_fz,csm_cz,csm_pz,csm_eog1,csm_eog2,csm_ntrig]]  
        rawEEG = rawData[:,[csm_fz,csm_cz,csm_pz]]              # rawData4 is rawData3 without trigger channel (index 5)
        rawEOG = rawData[:,[csm_eog1,csm_eog2]]
    
        # time period of scan
        fs = gtec['fs']
        sfreq = fs
        dt = 1/sfreq
        stop = dataRows/sfreq 
        Ts = np.arange(0,stop,dt) 
        Ts = Ts.reshape((len(Ts)),1)
    
        # Adaptive filter
        eegData = rawEEG
        eogData = rawEOG
        forgetF=0.995
        nKernel=5
        startSample=0
        p = False
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
        
        # %% EEG Signal Quality Metrics: first time domain metric using score 1 is based on the rms  
        fz_Q1 = general_amplitude(cleanData[:,0],2000,500)
        cz_Q1 = general_amplitude(cleanData[:,1],2000,500)
        pz_Q1 = general_amplitude(cleanData[:,2],2000,500) # (len(x.T)) is equal to the number of windows
        
        # %%  EEG Signal Quality Metrics: second time domain metric using score 3  : zero crossing rate
        fz_Q2 = amplitude_spectrum(cleanData[:,0],fs,2000,500)
        cz_Q2 = amplitude_spectrum(cleanData[:,1],fs,2000,500)
        pz_Q2 = amplitude_spectrum(cleanData[:,2],fs,2000,500)
        
        # %%  EEG Signal Quality Metrics: score 6: Kurtosis
        fz_Q3 = maximum_gradient(cleanData[:,0],2000,500)
        cz_Q3 = maximum_gradient(cleanData[:,1],2000,500)
        pz_Q3 = maximum_gradient(cleanData[:,2],2000,500)

        # Total score for fz channel
        fz_Quality = total_quality(fz_Q1[1], fz_Q2[1], fz_Q3[1])
        cz_Quality = total_quality(cz_Q1[1], cz_Q2[1], cz_Q3[1])
        pz_Quality = total_quality(pz_Q1[1], pz_Q2[1], pz_Q3[1])
        eegQuality = (np.array((fz_Quality,cz_Quality,pz_Quality)))

    return eegQuality

def sqf_gnd_avg(no_scans,no_channels,data,channel):
    metrics = np.array([])
    no_scans = no_scans
    no_chans = no_channels
    channel = channel
    folders = data
    for i in range(no_scans):
        eegQuality = sq_chan(folders[i])
        metrics = np.hstack([metrics, eegQuality])
    
    # place scores in respective channels
    metrics = metrics.reshape((len(metrics)),1)
    chan_qs = np.int_(metrics.reshape(no_scans,no_chans))
    
    # select channel cz for utilization
    scansQS = chan_qs[:,channel]
    # group the scans by either scan 1 or scan 2
    test_list = folders
    # using sorted() + groupby()
    # Initial Character Case Categorization
    util_func = lambda x: x[0:6]
    temp = sorted(test_list, key = util_func)
    runs = [list(ele) for i, ele in groupby(temp, util_func)]
    countScanTime = []
    for x in runs:
       p = (len(x))
       countScanTime.append(p)
       
    # first stage: group the quality scores 
    runsQS = []
    for number in countScanTime:
        runsQS.append(scansQS[:number])
        scansQS = scansQS[number:]
        
    # find mean of each element inside the list out
    avg_runsQS = []
    for i in range(len(runsQS)):
       avg_runsQS.append(np.mean(runsQS[i]))
    
    
    #%% extract the number of times group of runs belonging to a particular subject appear
    # automatically count the group of of scans belonging to participants
    xscans = []
    for i in range(len(runs)):
        xscans.append(runs[i][0])
    xscans = [x[:4] for x in xscans]
    count = Counter(map(tuple, xscans))
    count = pd.DataFrame.from_dict(count, orient='index').reset_index()
    countSubjRun = (count[0]).tolist()
    
    #%%
    # second stage: use the out_put to compare with the res file to scans per subject
    avg_runsQS_ = list(avg_runsQS)
    runs_ = runs
    subjQS = []
    subj = []
    for number in countSubjRun:
        subjQS.append(avg_runsQS_[:number])
        avg_runsQS_ = avg_runsQS_[number:]
        subj.append(runs_[:number])
        runs_ = runs_[number:]

    #%% maximum values
    
    # scores
    # find max values amongst grp_scores
    max_subjQS = []
    for i in range(len(subjQS)):
        if subjQS[0] == subjQS[1]:
            max_subjQS.append((subjQS[i][0]))
        else:
            max_subjQS.append(np.amax(subjQS[i]))
       
    # locate index of max values amongst grp_scores
    idx_ = []
    for i in range(len(subjQS)):
        idx_.append((np.where(subjQS[i] == max_subjQS[i]))[0])
        
    idx__ = []
    for i in range(len(subjQS)):
        idx__.append(idx_[i][0])
    
    # scans
    # use the idx of the maximum scores to extract the scans with maximum scores
    max_subj = []
    for i in range(len(idx__)):
        x = idx__[i]
        max_subj.append(subj[i][x])
    return max_subj


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