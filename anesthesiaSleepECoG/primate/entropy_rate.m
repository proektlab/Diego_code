
Skip to content
Pull requests
Issues
Marketplace
Explore
@diegodav
proektlab /
Diego_code
Public

Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights

    Settings

Diego_code/entropy/entropy_rate.m
@diegodav
diegodav Create entropy_rate.m
Latest commit 02142ae 5 days ago
History
1 contributor
23 lines (21 sloc) 547 Bytes
% Diego G. Davila
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This funciton calculates the entropy rate for a transition probability
% matrix. 
%
% Inputs: 
% 1. TRANS: Transition Probability Matrix
% 2. Mu: Steady State Probability Vector
%
% Outputs:
% 1. H: Entropy Rate of the Transition Probability Matrix
% 
function [H] = entropy_rate(TRANS, Mu)

H = 0; % initialize entropy as zero
for i = 1:length(Mu)
    H = H + nansum(Mu(i)*(TRANS(i,:).*log(TRANS(i,:)))); % sum each state's terms
end
H = -H; % flip the sign

end

    © 2022 GitHub, Inc.

    Terms
    Privacy
    Security
    Status
    Docs
    Contact GitHub
    Pricing
    API
    Training
    Blog
    About

Loading complete
