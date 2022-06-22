from fn_cfg import *
import params as cfg


filename = "No_Threshold.csv"
localPath_scans = "/Users/joshuaighalo/Documents/BrainNet/Projects/Workspace/results/laurel place/quality scans/"
file_dir = localPath_scans+filename
localPath_data = '/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/cleaned_dataset'
image_dest_dir = '/Users/joshuaighalo/Documents/BrainNet/Projects/Workspace/results/laurel place/images/grand_averages'
version = 1.0
disp = "Cz"
cutoff_LP = [False,False]


#   baseline scans for different dementia classes
quality_scans = pd.read_csv(file_dir)


base_ND = [x for x in (quality_scans['Baseline No Dementia'].to_numpy()) if str(x) != 'nan']
base_MILD = [x for x in (quality_scans['Baseline Mild Dementia'].to_numpy()) if str(x) != 'nan']
base_MODERATE = [x for x in (quality_scans['Baseline Moderate Dementia'].to_numpy()) if str(x) != 'nan']
base_SEVERE = [x for x in (quality_scans['Baseline Severe Dementia'].to_numpy()) if str(x) != 'nan']


four_ND = [x for x in (quality_scans['4-months No Dementia'].to_numpy()) if str(x) != 'nan']
four_MILD = [x for x in (quality_scans['4-months Mild Dementia'].to_numpy()) if str(x) != 'nan']
four_MODERATE = [x for x in (quality_scans['4-months Moderate Dementia'].to_numpy()) if str(x) != 'nan']
four_SEVERE = [x for x in (quality_scans['4-months Severe Dementia'].to_numpy()) if str(x) != 'nan']

eight_ND = [x for x in (quality_scans['8-months No Dementia'].to_numpy()) if str(x) != 'nan']
eight_MILD = [x for x in (quality_scans['8-months Mild Dementia'].to_numpy()) if str(x) != 'nan']
eight_MODERATE = [x for x in (quality_scans['8-months Moderate Dementia'].to_numpy()) if str(x) != 'nan']
eight_SEVERE = [x for x in (quality_scans['8-months Severe Dementia'].to_numpy()) if str(x) != 'nan']

print("ERPs of Cognitive Impairment Classes @ Baseline")
erp_baseND = averageERPs(device_version=version,scan_IDs=base_ND,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="No Dementia @ Baseline",img_name="No Dementia @ Baseline",destination_dir=image_dest_dir)
erp_baseMILD = averageERPs(device_version=version,scan_IDs=base_MILD,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Mild Dementia @ Baseline",img_name="Mild Dementia @ Baseline",destination_dir=image_dest_dir)
erp_baseMODERATE = averageERPs(device_version=version,scan_IDs=base_MODERATE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Moderate Dementia @ Baseline",img_name="Moderate Dementia @ Baseline",destination_dir=image_dest_dir)
erp_baseSEVERE = averageERPs(device_version=version,scan_IDs=base_SEVERE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Severe Dementia @ Baseline",img_name="Severe Dementia @ Baseline",destination_dir=image_dest_dir)


print("ERPs of Cognitive Impairment Classes @ 4-months")
erp_FourND = averageERPs(device_version=version,scan_IDs=four_ND,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="No Dementia @ 4-months",img_name="No Dementia @ 4-months",destination_dir=image_dest_dir)
erp_FourMILD = averageERPs(device_version=version,scan_IDs=four_MILD,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Mild Dementia @ 4-months",img_name="Mild Dementia @ 4-months",destination_dir=image_dest_dir)
erp_FourMODERATE = averageERPs(device_version=version,scan_IDs=four_MODERATE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Moderate Dementia @ 4-months",img_name="Moderate Dementia @ 4-months",destination_dir=image_dest_dir)
erp_FourSEVERE = averageERPs(device_version=version,scan_IDs=four_SEVERE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Severe Dementia @ 4-months",img_name="Severe Dementia @ 4-months",destination_dir=image_dest_dir)

print("ERPs of Cognitive Impairment Classes @ 8-months")
erp_EightND = averageERPs(device_version=version,scan_IDs=eight_ND,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="No Dementia @ 8-months",img_name="No Dementia @ 8-months",destination_dir=image_dest_dir)
erp_EightMILD = averageERPs(device_version=version,scan_IDs=eight_MILD,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Mild Dementia @ 8-months",img_name="Mild Dementia @ 8-months",destination_dir=image_dest_dir)
erp_EightMODERATE = averageERPs(device_version=version,scan_IDs=eight_MODERATE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Moderate Dementia @ 8-months",img_name="Moderate Dementia @ 8-months",destination_dir=image_dest_dir)
erp_EightSEVERE = averageERPs(device_version=version,scan_IDs=eight_SEVERE,dispIMG_Channel=disp,local_path=localPath_data,fs=cfg.fs,line=cfg.line,lowcut=cfg.highPass,highcut=cfg.lowPass,stimTrig=cfg.stimTrig,clip=cfg.clip,lowPassERP=cutoff_LP,label="Severe Dementia @ 8-months",img_name="Severe Dementia @ 8-months",destination_dir=image_dest_dir)


