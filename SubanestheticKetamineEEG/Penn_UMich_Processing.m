%% Penn UMich Joint Dataset Processing
%
% Diego G. Davila
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This script jointly processes EEG data from Penn and UMich Ketamine Studies
%
% Processing Steps:
% 1. Separate Data By Dose
% 2. Remove Non-Scalp Channels
% 3. Downsample Penn Data To Match UMich Data Sampling Rate
% 4. Apply a High-Pass Filter with Cutoff at 1Hz
% 5. Extract 1 minute of cleanest data to keep for further processing and
%    analysis
% 6. Remove and Interpolate Bad Channels
% 7. Perform ICA  
% 8. Remove line noise via CleanLine 
% 9. Automated Artifact Labeling & Rejection (Reguires this function: https://github.com/proektlab/Diego_code/blob/main/UsefulMiscFunctions/clean_comp.m)
% 10. Perform Average Re-Referencing 
%
% Required Packages & Functions: EEGLAB, Proekt Lab Spectral Analysis Directory 
addpath 'C:\Users\diego\Documents\MATLAB\eeglab2020_0';
eeglab % initialize eeglab
close all; % close eeglab window
addpath 'Y:\code\SpectralAnalysis';

%% Read in Raw Penn Data

% set up some useful lists
sl = [5 5 5 6 6 6 8 8 8 9 9 9 11 11 11 12 12 12 13 13 13]; % indexwise subject assignment
st = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3]; % indexwise dose assignment
subject_list = [5 6 8 9 11 12 13]; % list of all subjects

% read in the data
Penn_raw = {[]};
for i = 1:21
    Penn_raw{i} = pop_loadset(['sub' num2str(sl(i)) '_step' num2str(st(i)) '.set'], 'Y:\diego\ketamineConnectomics\EEGDATASETS\RAW');
end

% separate out the datasets by dosage
Penn_raw_baseline = Penn_raw(1:3:21); % no ketamine
Penn_raw_low_dose = Penn_raw(2:3:21); % low subanesthetic dose
Penn_raw_high_dose = Penn_raw(3:3:21); % high subanesthetic dose

%% Read in Raw UMich Data

UMich_raw = {[]};
for i = 1:15
    if (i<10)
        UMich_raw{i} = pop_loadset(['kk_0',num2str(i),'_Raw.set'], 'Y:\diego\ketamineConnectomics\Mashour_Data');
    else
        UMich_raw{i} = pop_loadset(['kk_',num2str(i),'_Raw.set'], 'Y:\diego\ketamineConnectomics\Mashour_Data');
    end
end

%% Extract UMich Subanesthetic Period

% define what the subanesthetic sample points are (Given by data sheet)
UMich_suban_periods = {[1624458 2816024], 
    [567281 1749724],
    [695379	1841172],
    [648882	1852448],
    [473066	1670080],
    [568583	1764477],
    [768602	1972072],
    [855071	2003940],
    [915496	2106663],
    [1015523 2222716],
    [516459	1714678],
    [784180	1989495],
    [1336040 2520605],
    [375161	1560230],
    [852600	2039694]};

UMich_suban_raw = UMich_raw; % prepare a chunked copy
for i = 1:length(UMich_suban_raw)
    bound1 = [0 UMich_suban_periods{i}(1)]; % lower time bound condition
    bound2 = [UMich_suban_periods{i}(2) length(UMich_suban_raw{i}.data(1,:))]; % upper time bound condition
    UMich_suban_raw{i} = eeg_eegrej(UMich_suban_raw{i}, bound2); % cut the top off first
    UMich_suban_raw{i} = eeg_eegrej(UMich_suban_raw{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

%% Extract UMich Baseline Period

% define what the baseline sample points are (Given by data sheet)
UMich_baseline_periods = {[597006	752534], 
    [364176 519705],
    [524021	679549],
    [242596	398124],
    [262594	418122],
    [410139	565667],
    [604638	760166],
    [626942	782470],
    [746031	901559],
    [777808	933337],
    [339124	497183],
    [536024	694084],
    [1113983 1272052],
    [163491	321553],
    [650994	809057]};

UMich_baseline_raw = UMich_raw; % prepare a chunked copy
for i = 1:length(UMich_baseline_raw)
    bound1 = [0 UMich_baseline_periods{i}(1)]; % lower time bound condition
    bound2 = [UMich_baseline_periods{i}(2) length(UMich_baseline_raw{i}.data(1,:))]; % upper time bound condition
    UMich_baseline_raw{i} = eeg_eegrej(UMich_baseline_raw{i}, bound2); % cut the top off first
    UMich_baseline_raw{i} = eeg_eegrej(UMich_baseline_raw{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

%% Remove non-scalp channels
% Both datasets used same channel mapping, so we can just define which
% channels to keep

% Scalp Channels 
scalp_chans = [15,22,9,26,23,18,16,10,3,2,33,27,24,19,11,4,124,123,34,28,20,12,5,118,117,116,40,35,29,13,6,112,111,110,109,41,36,30,7,106,105,104,103,46,47,42,37,31,80,87,93,98,102,51,52,53,54,55,79,86,92,97,58,59,60,61,78,85,91,65,66,67,62,77,84,90,70,71,72,76,83,75];

% Remove non-scalp channels from Penn data
for i = 1:length(Penn_raw_baseline)
    data = pop_select(Penn_raw_baseline{i}, 'channel', scalp_chans); % remove the non-scalp channels from baseline dataset
    Penn_raw_baseline{i} = data; % overwrite
    
    data = pop_select(Penn_raw_low_dose{i}, 'channel', scalp_chans); % remove the non-scalp channels from low dose dataset
    Penn_raw_low_dose{i} = data; % overwrite
    
    data = pop_select(Penn_raw_high_dose{i}, 'channel', scalp_chans); % remove the non-scalp channels from high dose dataset
    Penn_raw_high_dose{i} = data; % overwrite
end

% Remove non-scalp channels from UMich data
for i = 1:length(UMich_baseline_raw)
    data = pop_select(UMich_baseline_raw{i}, 'channel', scalp_chans); % remove the non-scalp channels from baseline dataset
    UMich_baseline_raw{i} = data; % overwrite

    data = pop_select(UMich_suban_raw{i}, 'channel', scalp_chans); % remove the non-scalp channels from ketamine dataset
    UMich_suban_raw{i} = data; % overwrite
end

%% Downsample Penn Data to Match UMich Sampling Rate
new_resolution = 500; % set the new sampling rate in Hz

for i = 1:length(Penn_raw_baseline)
    
    % baseline data
    Penn_raw_baseline{i} = pop_resample(Penn_raw_baseline{i}, new_resolution);
    % low dose data
    Penn_raw_low_dose{i} = pop_resample(Penn_raw_low_dose{i}, new_resolution);
    % high dose data
    Penn_raw_high_dose{i} = pop_resample(Penn_raw_high_dose{i}, new_resolution);
    
end

%% Filter Data
% Apply High-Pass Filtering

filt_cutoff = 1; % set low-end filter cutoff frequency in Hz

% Penn Data
for i = 1:length(Penn_raw_baseline)
    % baseline dose 
    data = pop_eegfiltnew(Penn_raw_baseline{i}, 'locutoff', filt_cutoff); % filter
    Penn_raw_baseline{i} = data; % overwrite
    
    % low dose
    data = pop_eegfiltnew(Penn_raw_low_dose{i}, 'locutoff', filt_cutoff); % filter
    Penn_raw_low_dose{i} = data; % overwrite
    
    % high dose
    data = pop_eegfiltnew(Penn_raw_high_dose{i}, 'locutoff', filt_cutoff); % filter
    Penn_raw_high_dose{i} = data; % overwrite
end

% UMich Data
for i = 1:length(UMich_baseline_raw)
    
    % baseline
    data = pop_eegfiltnew(UMich_baseline_raw{i}, 'locutoff', filt_cutoff); % filter
    UMich_baseline_raw{i} = data; % overwrite
    
    % ketamine
    data = pop_eegfiltnew(UMich_suban_raw{i}, 'locutoff', filt_cutoff); % filter
    UMich_suban_raw{i} = data; % overwrite
end

%% Select segments of time to use
% Extract 1 min of data. To get the corresponding indexes, we'll set the time intervals in seconds and
% multiply by the sampling rate as such: [start, stop] * sampling_rate

s_rate = 500; % sampling rate in Hz

Penn_baseline_times = {[15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate};

Penn_low_dose_times = {[15, 75] * s_rate,
    [15, 75] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate,
    [60, 120] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate};

Penn_high_dose_times = {[30, 90] * s_rate,
    [60, 120] * s_rate,
    [30, 90] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate};

UMich_baseline_times = {[75, 135] * s_rate,
    [15, 75] * s_rate,
    [120, 180] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [30, 90] * s_rate,
    [60, 120] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate,
    [60, 120] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate,
    [15, 75] * s_rate};

UMich_suban_times = {[30, 90] * s_rate,
    [30, 90] * s_rate,
    [240, 300] * s_rate,
    [30, 90] * s_rate,
    [330, 390] * s_rate,
    [120, 180] * s_rate,
    [120, 180] * s_rate,
    [30, 90] * s_rate,
    [30, 90] * s_rate,
    [60, 120] * s_rate,
    [30, 90] * s_rate,
    [15, 75] * s_rate,
    [60, 120] * s_rate,
    [60, 120] * s_rate,
    [60, 120] * s_rate};

%% Extract segments of time defined in previous step

% Penn Baseline
Penn_processed_baseline = Penn_raw_baseline; % prepare a copy
for i = 1:length(Penn_processed_baseline)
    bound1 = [0 Penn_baseline_times{i}(1)]; % lower time bound condition
    bound2 = [Penn_baseline_times{i}(2) length(Penn_processed_baseline{i}.data(1,:))]; % upper time bound condition
    Penn_processed_baseline{i} = eeg_eegrej(Penn_processed_baseline{i}, bound2); % cut the top off first
    Penn_processed_baseline{i} = eeg_eegrej(Penn_processed_baseline{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

% Penn Low Dose
Penn_processed_low_dose = Penn_raw_low_dose; % prepare a copy
for i = 1:length(Penn_processed_low_dose)
    bound1 = [0 Penn_low_dose_times{i}(1)]; % lower time bound condition
    bound2 = [Penn_low_dose_times{i}(2) length(Penn_processed_low_dose{i}.data(1,:))]; % upper time bound condition
    Penn_processed_low_dose{i} = eeg_eegrej(Penn_processed_low_dose{i}, bound2); % cut the top off first
    Penn_processed_low_dose{i} = eeg_eegrej(Penn_processed_low_dose{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

% Penn High Dose
Penn_processed_high_dose = Penn_raw_high_dose; % prepare a copy
for i = 1:length(Penn_processed_high_dose)
    bound1 = [0 Penn_high_dose_times{i}(1)]; % lower time bound condition
    bound2 = [Penn_high_dose_times{i}(2) length(Penn_processed_high_dose{i}.data(1,:))]; % upper time bound condition
    Penn_processed_high_dose{i} = eeg_eegrej(Penn_processed_high_dose{i}, bound2); % cut the top off first
    Penn_processed_high_dose{i} = eeg_eegrej(Penn_processed_high_dose{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

% UMich Baseline
UMich_baseline_processed = UMich_baseline_raw; % prepare a copy
for i = 1:length(UMich_baseline_processed)
    bound1 = [0 UMich_baseline_times{i}(1)]; % lower time bound condition
    bound2 = [UMich_baseline_times{i}(2) length(UMich_baseline_processed{i}.data(1,:))]; % upper time bound condition
    UMich_baseline_processed{i} = eeg_eegrej(UMich_baseline_processed{i}, bound2); % cut the top off first
    UMich_baseline_processed{i} = eeg_eegrej(UMich_baseline_processed{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

% UMich Subanesthetic Dose
UMich_suban_processed = UMich_suban_raw; % prepare a copy
for i = 1:length(UMich_suban_processed)
    bound1 = [0 UMich_suban_times{i}(1)]; % lower time bound condition
    bound2 = [UMich_suban_times{i}(2) length(UMich_suban_processed{i}.data(1,:))]; % upper time bound condition
    UMich_suban_processed{i} = eeg_eegrej(UMich_suban_processed{i}, bound2); % cut the top off first
    UMich_suban_processed{i} = eeg_eegrej(UMich_suban_processed{i}, bound1); % then, cut the bottom (this order prevents indexing errors)
end

%% Remove and interpolate bad channels
% remove channels that are std_thresh standard deviations from the mean

std_thresh = 3; % define threshold for rejection (i.e. reject channels this many standard deviations from the mean)

% Penn Baseline data
for i = 1:length(Penn_processed_baseline)
    data = Penn_processed_baseline{i}; % get the dataset
    [temp_set, bad_ind] = pop_rejchanspec(data, 'stdthresh', std_thresh); % get indexes of channels 2 standard deviations away from mean
    if isempty(bad_ind) == false % only interpolate if any bad channels were identified
        new_set = eeg_interp(data, bad_ind, 'spherical'); % spherical spline interpolation, using indexes derived just above
        Penn_processed_baseline{i} = new_set; % overwrite
    end
end

% Penn Low Dose data
for i = 1:length(Penn_processed_low_dose)
    data = Penn_processed_low_dose{i}; % get the dataset
    [temp_set, bad_ind] = pop_rejchanspec(data, 'stdthresh', std_thresh); % get indexes of channels 2 standard deviations away from mean
    if isempty(bad_ind) == false % only interpolate if any bad channels were identified
        new_set = eeg_interp(data, bad_ind, 'spherical'); % spherical spline interpolation, using indexes derived just above
        Penn_processed_low_dose{i} = new_set; % overwrite
    end
end

% Penn High Dose data
for i = 1:length(Penn_processed_high_dose)
    data = Penn_processed_high_dose{i}; % get the dataset
    [temp_set, bad_ind] = pop_rejchanspec(data, 'stdthresh', std_thresh); % get indexes of channels 2 standard deviations away from mean
    if isempty(bad_ind) == false % only interpolate if any bad channels were identified
        new_set = eeg_interp(data, bad_ind, 'spherical'); % spherical spline interpolation, using indexes derived just above
        Penn_processed_high_dose{i} = new_set; % overwrite
    end
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    data = UMich_baseline_processed{i}; % get the dataset
    [temp_set, bad_ind] = pop_rejchanspec(data, 'stdthresh', std_thresh); % get indexes of channels 2 standard deviations away from mean
    if isempty(bad_ind) == false % only interpolate if any bad channels were identified
        new_set = eeg_interp(data, bad_ind, 'spherical'); % spherical spline interpolation, using indexes derived just above
        UMich_baseline_processed{i} = new_set; % overwrite
    end
end

% UMich Subanesthetic Data
for i = 1:length(UMich_suban_processed)
    data = UMich_suban_processed{i}; % get the dataset
    [temp_set, bad_ind] = pop_rejchanspec(data, 'stdthresh', std_thresh); % get indexes of channels 2 standard deviations away from mean
    if isempty(bad_ind) == false % only interpolate if any bad channels were identified
        new_set = eeg_interp(data, bad_ind, 'spherical'); % spherical spline interpolation, using indexes derived just above
        UMich_suban_processed{i} = new_set; % overwrite
    end
end

%% Run ICA 
% Run ICA on all of the sessions individually so we can use components to
% remove line noise and reject artifacts

% Penn Baseline data
for i = 1:length(Penn_processed_baseline)
    Penn_processed_baseline{i} = pop_runica(Penn_processed_baseline{i}, 'icatype', 'runica');
end

% Penn Low Dose data
for i = 1:length(Penn_processed_low_dose)
    Penn_processed_low_dose{i} = pop_runica(Penn_processed_low_dose{i}, 'icatype', 'runica');
end

% Penn High Dose data
for i = 1:length(Penn_processed_high_dose)
    Penn_processed_high_dose{i} = pop_runica(Penn_processed_high_dose{i}, 'icatype', 'runica');
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    UMich_baseline_processed{i} = pop_runica(UMich_baseline_processed{i}, 'icatype', 'runica');
end

% UMich Subanesthetic Data
for i = 1:length(UMich_suban_processed)
    UMich_suban_processed{i} = pop_runica(UMich_suban_processed{i}, 'icatype', 'runica');
end

% CALCULATE ICA ACTIVATIONS 

% Penn Baseline data
for i = 1:length(Penn_processed_baseline)
    proc_dat = Penn_processed_baseline{i}; % for convenience, make a copy
    Penn_processed_baseline{i}.icaact = icaact(proc_dat.data, proc_dat.icaweights*proc_dat.icasphere, mean(proc_dat.data'));
end

% Penn Low Dose data
for i = 1:length(Penn_processed_low_dose)
    proc_dat = Penn_processed_low_dose{i}; % for convenience, make a copy
    Penn_processed_low_dose{i}.icaact = icaact(proc_dat.data, proc_dat.icaweights*proc_dat.icasphere, mean(proc_dat.data'));
end

% Penn High Dose data
for i = 1:length(Penn_processed_high_dose)
    proc_dat = Penn_processed_high_dose{i}; % for convenience, make a copy
    Penn_processed_high_dose{i}.icaact = icaact(proc_dat.data, proc_dat.icaweights*proc_dat.icasphere, mean(proc_dat.data'));
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    proc_dat = UMich_baseline_processed{i}; % for convenience, make a copy
    UMich_baseline_processed{i}.icaact = icaact(proc_dat.data, proc_dat.icaweights*proc_dat.icasphere, mean(proc_dat.data'));
end

% UMich Subanesthetic Data
for i = 1:length(UMich_suban_processed)
    proc_dat = UMich_suban_processed{i}; % for convenience, make a copy
    UMich_suban_processed{i}.icaact = icaact(proc_dat.data, proc_dat.icaweights*proc_dat.icasphere, mean(proc_dat.data'));
end

%% Remove Line Noise
% Use cleanline pluggin for EEGLAB to remove line noise around 60 Hz and
% 120 Hz and 180 Hz harmonics. 

bd_thresh = 2; % Set the bandwidth for noise removal in Hz
noise_freq = [60 120 180]; % set frequencies to center around

% Penn Baseline data
for i = 1:length(Penn_processed_baseline)
    Penn_processed_baseline{i} = pop_cleanline(Penn_processed_baseline{i}, 'LineFrequencies', noise_freq, 'Bandwidth', bd_thresh, 'ComputeSpectralPower', 'false', 'VerboseOutput', 'false');
end

% Penn Low Dose data
for i = 1:length(Penn_processed_low_dose)
    Penn_processed_low_dose{i} = pop_cleanline(Penn_processed_low_dose{i}, 'LineFrequencies', noise_freq, 'Bandwidth', bd_thresh, 'ComputeSpectralPower', 'false', 'VerboseOutput', 'false');
end

% Penn High Dose data
for i = 1:length(Penn_processed_high_dose)
    Penn_processed_high_dose{i} = pop_cleanline(Penn_processed_high_dose{i}, 'LineFrequencies', noise_freq, 'Bandwidth', bd_thresh, 'ComputeSpectralPower', 'false', 'VerboseOutput', 'false');
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    UMich_baseline_processed{i} = pop_cleanline(UMich_baseline_processed{i}, 'LineFrequencies', noise_freq, 'Bandwidth', bd_thresh, 'ComputeSpectralPower', 'false', 'VerboseOutput', 'false');
end

% UMich Subanesthetic Data
for i = 1:length(UMich_suban_processed)
    UMich_suban_processed{i} = pop_cleanline(UMich_suban_processed{i}, 'LineFrequencies', noise_freq, 'Bandwidth', bd_thresh, 'ComputeSpectralPower', 'false', 'VerboseOutput', 'false');
end

%% Automated Component Labeling via ICLabel & Rejection of Non-Brain (Artifact) Components

source_likelihood_threshold = 0.7; % set criteria for likelihood of component coming from brain source (the algorithm must be at least this percent confident the component is from a brain source)

% Penn Baseline Data
for i = 1:length(Penn_processed_baseline)
    Penn_processed_baseline{i} = clean_comp(Penn_processed_baseline{i}, source_likelihood_threshold);
end

% Penn Low Dose Data
for i = 1:length(Penn_processed_low_dose)
    Penn_processed_low_dose{i} = clean_comp(Penn_processed_low_dose{i}, source_likelihood_threshold);
end

% Penn High Dose Data
for i = 1:length(Penn_processed_high_dose)
    Penn_processed_high_dose{i} = clean_comp(Penn_processed_high_dose{i}, source_likelihood_threshold);
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    UMich_baseline_processed{i} = clean_comp(UMich_baseline_processed{i}, source_likelihood_threshold);
end

% UMich Suban Data
for i = 1:length(UMich_suban_processed)
    UMich_suban_processed{i} = clean_comp(UMich_suban_processed{i}, source_likelihood_threshold);
end

%% Perform Average Re-Referencing

% Penn Baseline data
for i = 1:length(Penn_processed_baseline)
    data = pop_reref(Penn_processed_baseline{i},[]); % re-reference to average
    Penn_processed_baseline{i} = data; % overwrite
end

% Penn Low Dose data
for i = 1:length(Penn_processed_low_dose)
    data = pop_reref(Penn_processed_low_dose{i},[]); % re-reference to average
    Penn_processed_low_dose{i} = data; % overwrite
end

% Penn High Dose data
for i = 1:length(Penn_processed_high_dose)
    data = pop_reref(Penn_processed_high_dose{i},[]); % re-reference to average
    Penn_processed_high_dose{i} = data; % overwrite
end

% UMich Baseline Data
for i = 1:length(UMich_baseline_processed)
    data = pop_reref(UMich_baseline_processed{i},[]); % re-reference to average
    UMich_baseline_processed{i} = data; % overwrite
end

% UMich Subanesthetic Data
for i = 1:length(UMich_suban_processed)
    data = pop_reref(UMich_suban_processed{i},[]); % re-reference to average
    UMich_suban_processed{i} = data; % overwrite
end

%% Save Processed Data

save('ProcessedEEGDatasets.mat', 'Penn_processed_baseline', 'Penn_processed_low_dose', 'Penn_processed_high_dose', 'UMich_baseline_processed', 'UMich_suban_processed', '-v7.3');



