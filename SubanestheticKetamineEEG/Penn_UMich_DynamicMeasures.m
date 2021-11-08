%% Penn UMich Dynamic Measure Derivation
%
% Diego G. Davila
% NGG
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This script derives the following measures:
%
%
addpath 'C:\Users\diego\Documents\MATLAB\eeglab2020_0';
eeglab % initialize eeglab
close all; % close eeglab window
addpath 'Y:\code\SpectralAnalysis';
addpath 'C:\Users\diego\Documents\MATLAB\fieldtrip-lite-20210311\fieldtrip-20210311';
%load('Y:\diego\ketamineConnectomics\EEGDATASETS\PROCESSED\ProcessedEEGDatasets.mat') % Load Processed Data
%
%% Convert to FieldTrip Format
% Convert Penn Data
Penn_baseline_ft = {[]};
Penn_low_dose_ft = {[]};
Penn_high_dose_ft = {[]};
for i = 1:length(Penn_processed_baseline)
    Penn_baseline_ft{i} = eeglab2fieldtrip(Penn_processed_baseline{i}, 'raw', 'none'); % convert to FieldTrip format and save 
    Penn_low_dose_ft{i} = eeglab2fieldtrip(Penn_processed_low_dose{i}, 'raw', 'none'); % convert to FieldTrip format and save 
    Penn_high_dose_ft{i} = eeglab2fieldtrip(Penn_processed_high_dose{i}, 'raw', 'none'); % convert to FieldTrip format and save 
end

% Convert UMich Data
UMich_baseline_ft = {[]};
UMich_suban_ft = {[]};
for i = 1:length(UMich_baseline_processed)
    UMich_baseline_ft{i} = eeglab2fieldtrip(UMich_baseline_processed{i}, 'raw', 'none'); % convert to FieldTrip format and save 
    UMich_suban_ft{i} = eeglab2fieldtrip(UMich_suban_processed{i}, 'raw', 'none'); % convert to FieldTrip format and save 
end

%% Segment processed data into time bins
% Use the trials field in the FieldTrip structure to store each segment

tau = 5; % set the length of the time bin in seconds
tau = tau*500; % convert to 500Hz sampling rate

% Define segment boundaries
segments = [1, tau+1, tau*2+1, tau*3+1, tau*4+1, tau*5+1, tau*6+1, tau*7+1, tau*8+1, tau*9+1, tau*10+1, tau*11+1; tau, tau*2, tau*3, tau*4, tau*5, tau*6, tau*7, tau*8, tau*9, tau*10, tau*11, tau*12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'; 

% Divide Penn data into segments
Penn_baseline_segmented = {[]};
Penn_low_dose_segmented = {[]};
Penn_high_dose_segmented = {[]};
for i = 1:length(Penn_baseline_ft)
    cfg = [];
    cgf.trl = segments;
    Penn_baseline_segmented{i} = ft_redefinetrial(cgf, Penn_baseline_ft{i});
    Penn_low_dose_segmented{i} = ft_redefinetrial(cgf, Penn_low_dose_ft{i});
    Penn_high_dose_segmented{i} = ft_redefinetrial(cgf, Penn_high_dose_ft{i});
end

% Divide UMich data into segments
UMich_baseline_segmented = {[]};
UMich_suban_segmented = {[]};
for i = 1:length(UMich_baseline_ft)
    cfg = [];
    cgf.trl = segments;
    UMich_baseline_segmented{i} = ft_redefinetrial(cgf, UMich_baseline_ft{i});
    UMich_suban_segmented{i} = ft_redefinetrial(cgf, UMich_suban_ft{i});
end
    
%% Extract alpha

% Penn alpha
Penn_baseline_alpha_segmented = {[]}; 
Penn_low_dose_alpha_segmented = {[]}; 
Penn_high_dose_alpha_segmented = {[]}; 
for i=1:length(Penn_baseline_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 3; % make it plus or minus 3 Hz
        cfg.foi = 10; % 10 Hz Alpha 
        cfg.trials = s;
        % baseline
        Penn_baseline_alpha_segmented{i}{s} = ft_freqanalysis(cfg, Penn_baseline_segmented{i});
        % low dose
        Penn_low_dose_alpha_segmented{i}{s} = ft_freqanalysis(cfg, Penn_low_dose_segmented{i});
        % high dose
        Penn_high_dose_alpha_segmented{i}{s} = ft_freqanalysis(cfg, Penn_high_dose_segmented{i});
    end
end

% UMich alpha
UMich_baseline_alpha_segmented = {[]}; 
UMich_suban_alpha_segmented = {[]}; 
for i=1:length(UMich_baseline_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 3; % make it plus or minus 3 Hz
        cfg.foi = 10; % 10 Hz Alpha 
        cfg.trials = s;
        % baseline
        UMich_baseline_alpha_segmented{i}{s} = ft_freqanalysis(cfg, UMich_baseline_segmented{i});
        % subanesthetic
        UMich_suban_alpha_segmented{i}{s} = ft_freqanalysis(cfg, UMich_suban_segmented{i});
    end
end
%% Calculate iCoh

% Penn iCoh
Penn_baseline_iCoh_segmented = {[]}; 
Penn_low_dose_iCoh_segmented = {[]}; 
Penn_high_dose_iCoh_segmented = {[]};
for i=1:length(Penn_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        
        % baseline
        data = ft_connectivityanalysis(cfg, Penn_baseline_alpha_segmented{i}{s});
        Penn_baseline_iCoh_segmented{i}(:,:,s) = data.cohspctrm;
        % low dose
        data = ft_connectivityanalysis(cfg, Penn_low_dose_alpha_segmented{i}{s});
        Penn_low_dose_iCoh_segmented{i}(:,:,s) = data.cohspctrm;
        % high dose
        data = ft_connectivityanalysis(cfg, Penn_high_dose_alpha_segmented{i}{s});
        Penn_high_dose_iCoh_segmented{i}(:,:,s) = data.cohspctrm;
    end
end

% UMich iCoh
UMich_baseline_iCoh_segmented = {[]}; 
UMich_suban_iCoh_segmented = {[]}; 
for i=1:length(UMich_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        
        % baseline
        data = ft_connectivityanalysis(cfg, UMich_baseline_alpha_segmented{i}{s});
        UMich_baseline_iCoh_segmented{i}(:,:,s) = data.cohspctrm;
        % subanesthetic
        data = ft_connectivityanalysis(cfg, UMich_suban_alpha_segmented{i}{s});
        UMich_suban_iCoh_segmented{i}(:,:,s) = data.cohspctrm;
    end
end

%% Calculate Phase-Slope Index (http://doc.ml.tu-berlin.de/causality/psi.pdf)

% The operations here are fairly computationally intensive, so probably
% best to use parallel computing

% parpool(14) % open a 14-CPU parallel computing pool

% Calculate the Cross-Spectral Desnity
% Penn CSD
Penn_baseline_csd_segmented = {[]}; 
Penn_low_dose_csd_segmented = {[]}; 
Penn_high_dose_csd_segmented = {[]};
parfor i=1:length(Penn_baseline_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'powandcsd';
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 3; % make bands plus or minus 3 Hz
        cfg.foirange = [1 48]; % 1 Hz - 48 Hz
        cfg.trials = s;
        
        % baseline
        data = ft_freqanalysis(cfg, Penn_baseline_segmented{i});
        Penn_baseline_csd_segmented{i}{s} = data;
        % low dose
        data = ft_freqanalysis(cfg, Penn_low_dose_segmented{i});
        Penn_low_dose_csd_segmented{i}{s} = data;
        % high dose
        data = ft_freqanalysis(cfg, Penn_high_dose_segmented{i});
        Penn_high_dose_csd_segmented{i}{s} = data;
    end
end

% UMich CSD
UMich_baseline_csd_segmented = {[]}; 
UMich_suban_csd_segmented = {[]}; 
parfor i=1:length(UMich_baseline_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'powandcsd';
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 3; % make bands plus or minus 3 Hz
        cfg.foirange = [1 48]; % 1 Hz - 48 Hz
        cfg.trials = s;
        
        % baseline
        data = ft_freqanalysis(cfg, UMich_baseline_segmented{i});
        UMich_baseline_csd_segmented{i}{s} = data;
        % subanesthetic
        data = ft_freqanalysis(cfg, UMich_suban_segmented{i});
        UMich_suban_csd_segmented{i}{s} = data;
    end
end

% Now, we can use the CSDs generated above to calculate PSI 

% Penn PSI
Penn_baseline_psi_segmented = {[]}; 
Penn_low_dose_psi_segmented = {[]}; 
Penn_high_dose_psi_segmented = {[]};
parfor i=1:length(Penn_baseline_csd_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'psi';
        cfg.bandwidth = 3;
        
        % baseline
        data = ft_connectivityanalysis(cfg, Penn_baseline_csd_segmented{i}{s});
        Penn_baseline_csd_segmented{i}{s} = []; % to conserve memory, wipe CSD
        Penn_baseline_psi_segmented{i}{s} = data;
        % low dose
        data = ft_connectivityanalysis(cfg, Penn_low_dose_csd_segmented{i}{s});
        Penn_low_dose_csd_segmented{i}{s} = []; % to conserve memory, wipe CSD
        Penn_low_dose_psi_segmented{i}{s} = data;
        % high dose
        data = ft_connectivityanalysis(cfg, Penn_high_dose_csd_segmented{i}{s});
        Penn_high_dose_csd_segmented{i}{s} = []; % to conserve memory, wipe CSD
        Penn_high_dose_psi_segmented{i}{s} = data;
    end
end

% UMich PSI
UMich_baseline_psi_segmented = {[]}; 
UMich_suban_psi_segmented = {[]}; 
parfor i=1:length(UMich_baseline_csd_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'psi';
        cfg.bandwidth = 3;
        
        % baseline
        data = ft_connectivityanalysis(cfg, UMich_baseline_csd_segmented{i}{s});
        UMich_baseline_csd_segmented{i}{s} = []; % to conserve memory, wipe CSD
        UMich_baseline_psi_segmented{i}{s} = data;
        % subanesthetic
        data = ft_connectivityanalysis(cfg, UMich_suban_csd_segmented{i}{s});
        UMich_suban_csd_segmented{i}{s} = []; % to conserve memory, wipe CSD
        UMich_suban_psi_segmented{i}{s} = data;
    end
end

% Ok, now we'll extract the alpha band PSI by averaging across the alpha
% frequency band (7-13 Hz) 
alpha_start = 36;
alpha_end = 66; 

% Extract Penn Alpha PSI
Penn_baseline_alpha_psi = {[]};
Penn_low_dose_alpha_psi = {[]};
Penn_high_dose_alpha_psi = {[]};
for i = 1:length(Penn_baseline_psi_segmented) % for every subject
    for s = 1:12 % for every time step 
        % baseline
        data = Penn_baseline_psi_segmented{i}{s}.psispctrm(:, [alpha_start:alpha_end]); % get the vectors representing the alpha band (rows = frequency, column = channel pair PSI values)
        data = mean(data, 2); % get the mean values across the bands
        Penn_baseline_alpha_psi{i}(:,:,s) = squareform(data); % save the data in matrix form
        % low dose
        data = Penn_low_dose_psi_segmented{i}{s}.psispctrm(:, [alpha_start:alpha_end]); % get the vectors representing the alpha band (rows = frequency, column = channel pair PSI values)
        data = mean(data, 2); % get the mean values across the bands
        Penn_low_dose_alpha_psi{i}(:,:,s) = squareform(data); % save the data in matrix form
        % high dose
        data = Penn_high_dose_psi_segmented{i}{s}.psispctrm(:, [alpha_start:alpha_end]); % get the vectors representing the alpha band (rows = frequency, column = channel pair PSI values)
        data = mean(data, 2); % get the mean values across the bands
        Penn_high_dose_alpha_psi{i}(:,:,s) = squareform(data); % save the data in matrix form
    end
end

% Extract UMich Alpha PSI
UMich_baseline_alpha_psi = {[]};
UMich_suban_alpha_psi = {[]};
for i = 1:length(UMich_baseline_psi_segmented) % for every subject
    for s = 1:12 % for every time step 
        % baseline
        data = UMich_baseline_psi_segmented{i}{s}.psispctrm(:, [alpha_start:alpha_end]); % get the vectors representing the alpha band (rows = frequency, column = channel pair PSI values)
        data = mean(data, 2); % get the mean values across the bands
        UMich_baseline_alpha_psi{i}(:,:,s) = squareform(data); % save the data in matrix form
        % subanesthetic
        data = UMich_suban_psi_segmented{i}{s}.psispctrm(:, [alpha_start:alpha_end]); % get the vectors representing the alpha band (rows = frequency, column = channel pair PSI values)
        data = mean(data, 2); % get the mean values across the bands
        UMich_suban_alpha_psi{i}(:,:,s) = squareform(data); % save the data in matrix form
    end
end

%% Calculate Degree

% Penn degree
Penn_baseline_iCoh_deg_segmented = {[]}; 
Penn_low_dose_iCoh_deg_segmented = {[]}; 
Penn_high_dose_iCoh_deg_segmented = {[]};
for i=1:length(Penn_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'degrees';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, Penn_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_baseline_iCoh_deg_segmented{i}(:,s) = data.degrees;
        %Penn_baseline_iCoh_deg_segmented{i}(:,s) = thresholded_weighted_degree(Penn_baseline_iCoh_segmented{i}(:,:,s), 0.1);
        % low dose
        data = ft_connectivityanalysis(cfg, Penn_low_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_low_dose_iCoh_deg_segmented{i}(:,s) = data.degrees;
        %Penn_low_dose_iCoh_deg_segmented{i}(:,s) = thresholded_weighted_degree(Penn_low_dose_iCoh_segmented{i}(:,:,s), 0.1);
        % high dose
        data = ft_connectivityanalysis(cfg, Penn_high_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_high_dose_iCoh_deg_segmented{i}(:,s) = data.degrees;
        %Penn_high_dose_iCoh_deg_segmented{i}(:,s) = thresholded_weighted_degree(Penn_high_dose_iCoh_segmented{i}(:,:,s), 0.1);
    end
end

% UMich degree
UMich_baseline_iCoh_deg_segmented = {[]}; 
UMich_suban_iCoh_deg_segmented = {[]}; 
for i=1:length(UMich_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'degrees';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, UMich_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_baseline_iCoh_deg_segmented{i}(:,s) = data.degrees;
        %UMich_baseline_iCoh_deg_segmented{i}(:,s) = thresholded_weighted_degree(UMich_baseline_iCoh_segmented{i}(:,:,s), 0.1); 
        % low dose
        data = ft_connectivityanalysis(cfg, UMich_suban_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_suban_iCoh_deg_segmented{i}(:,s) = data.degrees;
        %UMich_suban_iCoh_deg_segmented{i}(:,s) = thresholded_weighted_degree(UMich_suban_iCoh_segmented{i}(:,:,s), 0.1);
    end
end

%% Calculate Betweenness Centrality 

% Penn betweenness
Penn_baseline_iCoh_bet_segmented = {[]}; 
Penn_low_dose_iCoh_bet_segmented = {[]}; 
Penn_high_dose_iCoh_bet_segmented = {[]};
for i=1:length(Penn_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'betweenness';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, Penn_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_baseline_iCoh_bet_segmented{i}(:,s) = data.betweenness;
        % low dose
        data = ft_connectivityanalysis(cfg, Penn_low_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_low_dose_iCoh_bet_segmented{i}(:,s) = data.betweenness;
        % high dose
        data = ft_connectivityanalysis(cfg, Penn_high_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_high_dose_iCoh_bet_segmented{i}(:,s) = data.betweenness;
    end
end

% UMich betweenness
UMich_baseline_iCoh_bet_segmented = {[]}; 
UMich_suban_iCoh_bet_segmented = {[]}; 
for i=1:length(UMich_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'betweenness';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, UMich_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_baseline_iCoh_bet_segmented{i}(:,s) = data.betweenness;
        % low dose
        data = ft_connectivityanalysis(cfg, UMich_suban_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_suban_iCoh_bet_segmented{i}(:,s) = data.betweenness;
    end
end

%% Calculate Clustering Coefficient 

% Penn degree
Penn_baseline_iCoh_cc_segmented = {[]}; 
Penn_low_dose_iCoh_cc_segmented = {[]}; 
Penn_high_dose_iCoh_cc_segmented = {[]};
for i=1:length(Penn_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'clustering_coef';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, Penn_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_baseline_iCoh_cc_segmented{i}(:,s) = data.clustering_coef;
        % low dose
        data = ft_connectivityanalysis(cfg, Penn_low_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_low_dose_iCoh_cc_segmented{i}(:,s) = data.clustering_coef;
        % high dose
        data = ft_connectivityanalysis(cfg, Penn_high_dose_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        Penn_high_dose_iCoh_cc_segmented{i}(:,s) = data.clustering_coef;
    end
end

% UMich degree
UMich_baseline_iCoh_cc_segmented = {[]}; 
UMich_suban_iCoh_cc_segmented = {[]}; 
for i=1:length(UMich_baseline_alpha_segmented)
    for s = 1:12
        cfg = [];
        cfg.method = 'coh'; % coherency
        cfg.complex = 'absimag'; % only get imaginary part of coherence spectrum, limit effect of volume conduction as per Nolte et al., 2004
        cfg2           = [];
        cfg2.method    = 'clustering_coef';
        cfg2.parameter = 'cohspctrm';
        cfg2.threshold = .1;
        % baseline
        data = ft_connectivityanalysis(cfg, UMich_baseline_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_baseline_iCoh_cc_segmented{i}(:,s) = data.clustering_coef;
        % low dose
        data = ft_connectivityanalysis(cfg, UMich_suban_alpha_segmented{i}{s});
        data = ft_networkanalysis(cfg2, data);
        UMich_suban_iCoh_cc_segmented{i}(:,s) = data.clustering_coef;
    end
end


%% iCoh Coeff. of Var

% Penn CV
Penn_baseline_iCoh_CV = {[]};
Penn_low_dose_iCoh_CV = {[]};
Penn_high_dose_iCoh_CV = {[]};
for i = 1:length(Penn_baseline_iCoh_segmented) % for every subject
    for m = 1:size(Penn_baseline_iCoh_segmented{i},1) % for every channel
        for n = 1:size(Penn_baseline_iCoh_segmented{i},2) % for every channel
            % baseline
            Penn_baseline_iCoh_CV{i}(m, n) = std(Penn_baseline_iCoh_segmented{i}(m, n, :))/mean(Penn_baseline_iCoh_segmented{i}(m, n, :)); % Calculate CV for that channel pair across time
            % low dose
            Penn_low_dose_iCoh_CV{i}(m, n) = std(Penn_low_dose_iCoh_segmented{i}(m, n, :))/mean(Penn_low_dose_iCoh_segmented{i}(m, n, :)); % Calculate CV for that channel pair across time
            % high dose
            Penn_high_dose_iCoh_CV{i}(m, n) = std(Penn_high_dose_iCoh_segmented{i}(m, n, :))/mean(Penn_high_dose_iCoh_segmented{i}(m, n, :)); % Calculate CV for that channel pair across time
        end
    end
end

% UMich CV
UMich_baseline_iCoh_CV = {[]};
UMich_suban_iCoh_CV = {[]};
for i = 1:length(UMich_baseline_iCoh_segmented) % for every subject
    for m = 1:size(UMich_baseline_iCoh_segmented{i},1) % for every channel
        for n = 1:size(UMich_baseline_iCoh_segmented{i},2) % for every channel
            % baseline
            UMich_baseline_iCoh_CV{i}(m, n) = std(UMich_baseline_iCoh_segmented{i}(m, n, :))/mean(UMich_baseline_iCoh_segmented{i}(m, n, :)); % Calculate CV for that channel pair across time
            % subanesthetic
            UMich_suban_iCoh_CV{i}(m, n) = std(UMich_suban_iCoh_segmented{i}(m, n, :))/mean(UMich_suban_iCoh_segmented{i}(m, n, :)); % Calculate CV for that channel pair across time
        end
    end
end

%% Degree Coeff. of Var

% Penn CV
Penn_baseline_iCoh_deg_CV = {[]};
Penn_low_dose_iCoh_deg_CV = {[]};
Penn_high_dose_iCoh_deg_CV = {[]};
for i = 1:length(Penn_baseline_iCoh_deg_segmented) % for every subject
    for m = 1:size(Penn_baseline_iCoh_deg_segmented{i},1) % for every channel
        % baseline
        Penn_baseline_iCoh_deg_CV{i}(m) = std(Penn_baseline_iCoh_deg_segmented{i}(m, :))/mean(Penn_baseline_iCoh_deg_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % low dose
        Penn_low_dose_iCoh_deg_CV{i}(m) = std(Penn_low_dose_iCoh_deg_segmented{i}(m, :))/mean(Penn_low_dose_iCoh_deg_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % high dose
        Penn_high_dose_iCoh_deg_CV{i}(m) = std(Penn_high_dose_iCoh_deg_segmented{i}(m, :))/mean(Penn_high_dose_iCoh_deg_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

% UMich CV
UMich_baseline_iCoh_deg_CV = {[]};
UMich_suban_iCoh_deg_CV = {[]};
for i = 1:length(UMich_baseline_iCoh_segmented) % for every subject
    for m = 1:size(UMich_baseline_iCoh_segmented{i},1) % for every channel
        % baseline
        UMich_baseline_iCoh_deg_CV{i}(m) = std(UMich_baseline_iCoh_deg_segmented{i}(m, :))/mean(UMich_baseline_iCoh_deg_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % subanesthetic
        UMich_suban_iCoh_deg_CV{i}(m) = std(UMich_suban_iCoh_deg_segmented{i}(m, :))/mean(UMich_suban_iCoh_deg_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

%% Betweenness Coeff. of Var

% Penn CV
Penn_baseline_iCoh_bet_CV = {[]};
Penn_low_dose_iCoh_bet_CV = {[]};
Penn_high_dose_iCoh_bet_CV = {[]};
for i = 1:length(Penn_baseline_iCoh_bet_segmented) % for every subject
    for m = 1:size(Penn_baseline_iCoh_bet_segmented{i},1) % for every channel
        % baseline
        Penn_baseline_iCoh_bet_CV{i}(m) = std(Penn_baseline_iCoh_bet_segmented{i}(m, :))/mean(Penn_baseline_iCoh_bet_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % low dose
        Penn_low_dose_iCoh_bet_CV{i}(m) = std(Penn_low_dose_iCoh_bet_segmented{i}(m, :))/mean(Penn_low_dose_iCoh_bet_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % high dose
        Penn_high_dose_iCoh_bet_CV{i}(m) = std(Penn_high_dose_iCoh_bet_segmented{i}(m, :))/mean(Penn_high_dose_iCoh_bet_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

% UMich CV
UMich_baseline_iCoh_bet_CV = {[]};
UMich_suban_iCoh_bet_CV = {[]};
for i = 1:length(UMich_baseline_iCoh_bet_segmented) % for every subject
    for m = 1:size(UMich_baseline_iCoh_bet_segmented{i},1) % for every channel
        % baseline
        UMich_baseline_iCoh_bet_CV{i}(m) = std(UMich_baseline_iCoh_bet_segmented{i}(m, :))/mean(UMich_baseline_iCoh_bet_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % subanesthetic
        UMich_suban_iCoh_bet_CV{i}(m) = std(UMich_suban_iCoh_bet_segmented{i}(m, :))/mean(UMich_suban_iCoh_bet_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

%% Clustering Coefficient Coeff. of Var

% Penn CV
Penn_baseline_iCoh_cc_CV = {[]};
Penn_low_dose_iCoh_cc_CV = {[]};
Penn_high_dose_iCoh_cc_CV = {[]};
for i = 1:length(Penn_baseline_iCoh_cc_segmented) % for every subject
    for m = 1:size(Penn_baseline_iCoh_cc_segmented{i},1) % for every channel
        % baseline
        Penn_baseline_iCoh_cc_CV{i}(m) = std(Penn_baseline_iCoh_cc_segmented{i}(m, :))/mean(Penn_baseline_iCoh_cc_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % low dose
        Penn_low_dose_iCoh_cc_CV{i}(m) = std(Penn_low_dose_iCoh_cc_segmented{i}(m, :))/mean(Penn_low_dose_iCoh_cc_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % high dose
        Penn_high_dose_iCoh_cc_CV{i}(m) = std(Penn_high_dose_iCoh_cc_segmented{i}(m, :))/mean(Penn_high_dose_iCoh_cc_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

% UMich CV
UMich_baseline_iCoh_cc_CV = {[]};
UMich_suban_iCoh_cc_CV = {[]};
for i = 1:length(UMich_baseline_iCoh_cc_segmented) % for every subject
    for m = 1:size(UMich_baseline_iCoh_cc_segmented{i},1) % for every channel
        % baseline
        UMich_baseline_iCoh_cc_CV{i}(m) = std(UMich_baseline_iCoh_cc_segmented{i}(m, :))/mean(UMich_baseline_iCoh_cc_segmented{i}(m, :)); % Calculate CV for that channel pair across time
        % subanesthetic
        UMich_suban_iCoh_cc_CV{i}(m) = std(UMich_suban_iCoh_cc_segmented{i}(m, :))/mean(UMich_suban_iCoh_cc_segmented{i}(m, :)); % Calculate CV for that channel pair across time
    end
end

%% Merge Penn and UMich Data

% merge segmented iCoh
all_baseline_iCoh_segmented = {Penn_baseline_iCoh_segmented{:} UMich_baseline_iCoh_segmented{:}};
all_suban_iCoh_segmented = {Penn_low_dose_iCoh_segmented{:} UMich_suban_iCoh_segmented{:}};
    
% merge CV matrices
all_baseline_iCoh_CV = {Penn_baseline_iCoh_CV{:} UMich_baseline_iCoh_CV{:}};
all_suban_iCoh_CV = {Penn_low_dose_iCoh_CV{:} UMich_suban_iCoh_CV{:}};
    
% merge degree matrices
all_baseline_iCoh_deg_segmented = {Penn_baseline_iCoh_deg_segmented{:} UMich_baseline_iCoh_deg_segmented{:}};
all_suban_iCoh_deg_segmented = {Penn_low_dose_iCoh_deg_segmented{:} UMich_suban_iCoh_deg_segmented{:}};
    
% merge degree cv matrices
all_baseline_iCoh_deg_CV = {Penn_baseline_iCoh_deg_CV{:} UMich_baseline_iCoh_deg_CV{:}};
all_suban_iCoh_deg_CV = {Penn_low_dose_iCoh_deg_CV{:} UMich_suban_iCoh_deg_CV{:}};

% merge segmented alpha PSI
all_baseline_alpha_psi = {Penn_baseline_alpha_psi{:} UMich_baseline_alpha_psi{:}};
all_suban_alpha_psi = {Penn_low_dose_alpha_psi{:} UMich_suban_alpha_psi{:}};

% merge CV matrices PENN HIGH KET
all_baseline_iCoh_CV = {Penn_baseline_iCoh_CV{:}};
all_suban_iCoh_CV = {Penn_high_dose_iCoh_CV{:}};

% merge degree cv matrices PENN HIGH KET
all_baseline_iCoh_deg_CV = {Penn_baseline_iCoh_deg_CV{:}};
all_suban_iCoh_deg_CV = {Penn_high_dose_iCoh_deg_CV{:}};

% merge betweenness cv matrices PENN HIGH KET
all_baseline_iCoh_bet_CV = {Penn_baseline_iCoh_bet_CV{:}};
all_suban_iCoh_bet_CV = {Penn_high_dose_iCoh_bet_CV{:}};

% merge clustering coefficient cv matrices PENN HIGH KET
all_baseline_iCoh_cc_CV = {Penn_baseline_iCoh_cc_CV{:}};
all_suban_iCoh_cc_CV = {Penn_high_dose_iCoh_cc_CV{:}};
    

%% Save 

save('Dynamic_Measures.mat', 'Penn_baseline_iCoh_segmented', 'Penn_low_dose_iCoh_segmented', 'Penn_high_dose_iCoh_segmented', 'UMich_baseline_iCoh_segmented', 'UMich_suban_iCoh_segmented', 'Penn_baseline_alpha_psi', 'Penn_low_dose_alpha_psi', 'Penn_high_dose_alpha_psi', 'UMich_baseline_alpha_psi', 'UMich_suban_alpha_psi', '-v7.3')

