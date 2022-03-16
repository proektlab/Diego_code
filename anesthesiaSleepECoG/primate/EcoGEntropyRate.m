function output = EcoGEntropyRate(data, fs, num_clust, segment_len)
% Diego G. Davila
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This function conducts an entropy rate analysis on markov models fitted
% to segments of ECoG data. 
% 
% INPUT(s):
% 1. data: time by channel matrix containing ECoG recordings. 
% 2. fs: sampling rate of data.
% 3. num_clust: numbe rof clusters used for k-means estimate
% 4. segment_len: how long of a period, in seconds, to estimate markov models over
%
% OUTPUT(s):
% 1. output: struct containing output measures
%
[idx,C,~] = kmeans(data, num_clust, 'Distance', 'cosine', 'MaxIter', 1000); % use cosine similarity as clutering criteria, as this ignores polarity

tau = segment_len; % in seconds, how long each period 
tau = tau * fs; % convert to sample points 

% get the number of segments this tau will give us 
segments = size(data, 1)/tau;
if (mod(segments, 1) > 0)
    segments = fix(segments);
end


markovs = {[]}; % to save markov models over time
entropy_vector = zeros(1,segments);

% split the chain into segments
idx_split = {[]};
for i = 1:(size(data, 1)/tau)
    idx_split{i} = idx(tau*i-tau + 1:tau*i);
end

trace_vector = zeros(1,segments);
% generate markov models and calculate entropy rate
for j = 1:segments % for every division point in the data
    markov_model = struct([]);
    markov_model(1).chain = idx_split{j};
    markov_model.states = C;

    % generate Markov model and populate the Markov object
    [TRANS, EMIS] = hmmestimate(idx_split{j}, idx_split{j}); 
    [V,D] = eig(TRANS');
    markov_model.eigenvalues = D;
    markov_model.eigenvectors = V;
    markov_model.transition = TRANS;
    markov_model.emissions = EMIS;
    markov_model.trace = trace(TRANS); 

    trace_vector(j) = markov_model.trace;
    % Get steady state probabilities
    x = markov_model.eigenvectors(:,1)';
    markov_model.steady_state = x(:,1)./sum(x(:,1));
    
    % calculate entropy rate
    markov_model.entropy_rate = entropy_rate(markov_model.transition, markov_model.steady_state);
    entropy_vector(j) = markov_model.entropy_rate;
    markovs{j} = markov_model;
end

output = struct([]);
output(1).markovs = markovs;
output.entropy_vector = entropy_vector;
output.trace_vector = trace_vector;

% PLOT
figure; plot(1:(size(data, 1)/tau), entropy_vector'); xlabel('Time Segments'); ylabel('Entropy Rate'); title('Entropy Rate Over Time');


end
