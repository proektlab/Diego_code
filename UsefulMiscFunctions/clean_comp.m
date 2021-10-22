function EEG = clean_comp(EEG, brain_source_percentage)
% This function takes as input an EEGLAB dataset (inEEG) and labels
% Independent Components via the ICLabel pluggin, and returns the same dataset with non-brain source components removed. 
%
% By: Diego G. Davila 
%     Proekt Lab 
%     University of Pennsylvania School of Medicine
%     10/06/2021
% 
% INPUTS: 
%     1. EEG: An EEGLAB Dataset that has had ICA Component labeling via
%     ICLabel.
%     2. brain_source_percentage: the minimum likelihood that a component
%     comes form a brain source.
% 
% OUTPUTS:
%     1. EEG: The cleaned EEGLAB dataset with artifact components rejected.
%
EEG = iclabel(EEG); % apply independent component labeling via ICLabel EEGLAB pluggin
brain_source_idx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= brain_source_percentage); % get the indexes of components likely to come from brain source
EEG = pop_subcomp(EEG, brain_source_idx, 0, 1); % keep components that are likely coming from brain source, and reject the rest
end
