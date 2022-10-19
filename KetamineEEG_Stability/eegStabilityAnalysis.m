function output = eegStabilityAnalysis(data, srate, subname, win)
%
% MODIFIED FROM: https://github.com/proektlab/stereoEEG/blob/main/seegStability.m
%
% INPUT(s):
% 1. data: time x channel EEG data matrix
% 2. srate: sampling rate in Hz
% 3. subname: name of the output files, recommend to use subject number (without the .mat extension; i.e. 'subject1')
% 4. win: window size for VAR fitting
%
% OUTPUT(s):
% 1. output: struct containing:
%    coefficient matrices
%    noise matrices
%    coeff. matrix eigenvalues
%    histograms of |lambda|
%    coeff. matrix eigenvectors
%    subject ID

% parameters
lpc = []; % filter cutoffs (Hz; leave empty [] for no filter)
car = 0; % common-average re-reference? (0 or 1)
twn = win; % time window size (s)

%lbn = 0.05:0.05:1.05; % eigenvalue bins

lbn = 0.05:0.005:1.02; % eigenvalue bins

lth = 0.5; % eigenvalue thresholds (can list more than one)
smt = [2 10]; % smoothing factor [crit index bins, time bins] (leave empty [] for no smoothing)

% set some variables
cdat = data;
fs = srate;
N = size(cdat, 2);

% time bins
is = 1:round(twn*fs):size(cdat,1); % window start times (samples)
Nwin = length(is)-1;
tcenter = is(1:end-1)/fs+twn/2; % s


% prep output
output = struct([]);
output(1).coefficient_matrices = zeros(N, N, Nwin);
output.noise_matrices = zeros(N, N, Nwin);
output.lambdas = zeros(N, Nwin);
output.lambda_histograms = []; 
output.subject = subname;
output.eigenvectors = zeros(size(cdat, 2), size(cdat, 2), Nwin);


% multivariate AR model fit and eigendecomposition
lambhist = zeros(length(lbn)-1,Nwin);
for ii = 1:Nwin
    [~,A,C,~,~,~] = arfit(cdat(is(ii):is(ii+1)-1,:),1,1); % first order model
    [~,~,~,~,~,lambda, vects] = armodeFAST_modified(A,C);
    lambhist(:,ii) = histcounts(abs(lambda),lbn)';
    
    % Save AR Model internals
    output.coefficient_matrices(:,:,ii) = A;
    output.noise_matrices(:,:,ii) = C; 
    output.lambdas(:,ii) = lambda;
    output.eigenvectors(:,:,ii) = vects;
end
if ~isempty(smt), lambhist = smooth2a(lambhist,smt(1),smt(2)); end
output.lambda_histograms = lambhist; % save smoothes lambda histograms


end

