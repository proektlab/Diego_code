%% Ketamine EEG Variance of Change Analysis
%
% Diego G. Davila
% NGG
% Proekt Lab
% University of Pennsylvania School of Medicine
%
% This script contrains an analysis of the variance in change in
% connectivity under 0.4 micrograms/mL of Ketamine
% 
% This is done for a measure of functional connectivity (Imaginary
% Coherence) and effective connectivity (Phase Slope Index)
% 
%% Generate Vectors of Connectivity Flucutation 

% iCoh

include_michigan = false; % boolean - determine wether to include UMich Dataset (increases baseline and low dose sample size by 15)

% Baseline
% Generate fluctuations over time segments for every subject at baseline
if (include_michigan)
    BaselineFluctuations_icoh=zeros(numel(all_baseline_iCoh_segmented), 55);
    for j=1:numel(all_baseline_iCoh_segmented)
        s1 = all_baseline_iCoh_segmented{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            BaselineFluctuations_icoh(j,i) = norm(s1(:,:,i+1) - s1(:,:,i), 'fro');
            normvector(i+1)=norm(s1(:,:,i+1) ,'fro');
        end
        BaselineFluctuations_icoh(j,:)=BaselineFluctuations_icoh(j,:)/mean(normvector);
    end
else
    BaselineFluctuations_icoh=zeros(numel(Penn_baseline_iCoh_segmented), 55);
    for j=1:numel(Penn_baseline_iCoh_segmented)
        s1 = Penn_baseline_iCoh_segmented{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            BaselineFluctuations_icoh(j,i) = norm(s1(:,:,i+1) - s1(:,:,i), 'fro');
            normvector(i+1)=norm(s1(:,:,i+1) ,'fro');
        end
        BaselineFluctuations_icoh(j,:)=BaselineFluctuations_icoh(j,:)/mean(normvector);
    end
end
% Low Dose
if (include_michigan)
    % Generate fluctuations over time segments for every subject under low ketamine
    lowKFluctuations_icoh = zeros(numel(all_suban_iCoh_segmented), 55);
    for j = 1:numel(all_suban_iCoh_segmented)
        s1 = all_suban_iCoh_segmented{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            lowKFluctuations_icoh(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
            normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
        end
        lowKFluctuations_icoh(j,:) = lowKFluctuations_icoh(j,:)/mean(normvector);
    end
else
    % Generate fluctuations over time segments for every subject under low ketamine
    lowKFluctuations_icoh = zeros(numel(Penn_low_dose_iCoh_segmented), 55);
    for j = 1:numel(Penn_low_dose_iCoh_segmented)
        s1 = Penn_low_dose_iCoh_segmented{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            lowKFluctuations_icoh(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
            normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
        end
        lowKFluctuations_icoh(j,:) = lowKFluctuations_icoh(j,:)/mean(normvector);
    end
end
% High Dose
% Generate fluctuations over time segments for every subject under ketamine
highKFluctuations_icoh = zeros(numel(Penn_high_dose_iCoh_segmented), 55);
for j = 1:numel(Penn_high_dose_iCoh_segmented)
    s1 = Penn_high_dose_iCoh_segmented{j};
    normvector = zeros(56, 1);
    normvector(1) = norm(s1(:,:,1) ,'fro');
    for i = 1:size(s1, 3) - 1
        highKFluctuations_icoh(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
        normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
    end
    highKFluctuations_icoh(j,:) = highKFluctuations_icoh(j,:)/mean(normvector);
end

% PSI

% Baseline
if (include_michigan)
    % Generate fluctuations over time segments for every subject at baseline
    BaselineFluctuations_psi=zeros(numel(all_baseline_alpha_psi), 55);
    for j=1:numel(all_baseline_alpha_psi)
        s1 = all_baseline_alpha_psi{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            BaselineFluctuations_psi(j,i) = norm(s1(:,:,i+1) - s1(:,:,i), 'fro');
            normvector(i+1)=norm(s1(:,:,i+1) ,'fro');
        end
        BaselineFluctuations_psi(j,:)=BaselineFluctuations_psi(j,:)/mean(normvector);
    end
else
    % Generate fluctuations over time segments for every subject at baseline
    BaselineFluctuations_psi=zeros(numel(Penn_baseline_alpha_psi), 55);
    for j=1:numel(Penn_baseline_alpha_psi)
        s1 = Penn_baseline_alpha_psi{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            BaselineFluctuations_psi(j,i) = norm(s1(:,:,i+1) - s1(:,:,i), 'fro');
            normvector(i+1)=norm(s1(:,:,i+1) ,'fro');
        end
        BaselineFluctuations_psi(j,:)=BaselineFluctuations_psi(j,:)/mean(normvector);
    end
end
% Low Dose
if (include_michigan)
    lowKFluctuations_psi = zeros(numel(all_suban_alpha_psi), 55);
    for j = 1:numel(all_suban_alpha_psi)
        s1 = all_suban_alpha_psi{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            lowKFluctuations_psi(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
            normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
        end
        lowKFluctuations_psi(j,:) = lowKFluctuations_psi(j,:)/mean(normvector);
    end
else
% Generate fluctuations over time segments for every subject under low ketamine
    lowKFluctuations_psi = zeros(numel(Penn_low_dose_alpha_psi), 55);
    for j = 1:numel(Penn_low_dose_alpha_psi)
        s1 = Penn_low_dose_alpha_psi{j};
        normvector = zeros(56, 1);
        normvector(1) = norm(s1(:,:,1) ,'fro');
        for i = 1:size(s1, 3) - 1
            lowKFluctuations_psi(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
            normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
        end
        lowKFluctuations_psi(j,:) = lowKFluctuations_psi(j,:)/mean(normvector);
    end
end
% High Dose
% Generate fluctuations over time segments for every subject under ketamine
highKFluctuations_psi = zeros(numel(Penn_high_dose_alpha_psi), 55);
for j = 1:numel(Penn_high_dose_alpha_psi)
    s1 = Penn_high_dose_alpha_psi{j};
    normvector = zeros(56, 1);
    normvector(1) = norm(s1(:,:,1) ,'fro');
    for i = 1:size(s1, 3) - 1
        highKFluctuations_psi(j,i) = norm(s1(:,:,i+1)-s1(:,:,i), 'fro');
        normvector(i+1) = norm(s1(:,:,i+1) ,'fro');
    end
    highKFluctuations_psi(j,:) = highKFluctuations_psi(j,:)/mean(normvector);
end

%% Generate Bootstrapped Estimates of Variance

num_boots = 100000; % define number of bootstrap runs

% PSI
baseline_var_boot_psi = bootstrp(num_boots, @var, BaselineFluctuations_psi(:));
low_dose_var_boot_psi = bootstrp(num_boots, @var, lowKFluctuations_psi(:));
high_dose_var_boot_psi = bootstrp(num_boots, @var, highKFluctuations_psi(:));
% iCoh
baseline_var_boot_icoh = bootstrp(num_boots, @var, BaselineFluctuations_icoh(:));
low_dose_var_boot_icoh = bootstrp(num_boots, @var, lowKFluctuations_icoh(:));
high_dose_var_boot_icoh = bootstrp(num_boots, @var, highKFluctuations_icoh(:));

%% Generate Plots 
% Generate plots demonstrating distributions of bootstrapped estimates of the mean and
% variance of fluctuations across Baseline and Ketamine conditions 

% Phase Slope Index
figure;
ksdensity(bootstrp(num_boots, @mean, BaselineFluctuations_psi(:)), 'Function', 'pdf');
hold on;
ksdensity(bootstrp(num_boots, @mean, lowKFluctuations_psi(:)), 'Function', 'pdf');
hold on;
ksdensity(bootstrp(num_boots, @mean, highKFluctuations_psi(:)), 'Function', 'pdf');
title(['Boostrapped estimates of the Mean of Fluctuations PSI | ' num2str(num_boots) ' Runs'])
xlabel('Mean Estimates')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')

figure;
ksdensity(baseline_var_boot_psi, 'Function', 'pdf');
hold on;
ksdensity(low_dose_var_boot_psi, 'Function', 'pdf');
hold on;
ksdensity(high_dose_var_boot_psi, 'Function', 'pdf');
title(['Boostrapped estimates of the Variance of Fluctuations PSI | ' num2str(num_boots) ' Runs'])
xlabel('Variance Estimates')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')

% Imaginary Coherence
figure;
ksdensity(bootstrp(num_boots, @mean, BaselineFluctuations_icoh(:)), 'Function', 'pdf');
hold on;
ksdensity(bootstrp(num_boots, @mean, lowKFluctuations_icoh(:)), 'Function', 'pdf');
hold on;
ksdensity(bootstrp(num_boots, @mean, highKFluctuations_icoh(:)), 'Function', 'pdf');
title(['Boostrapped estimates of the Mean of Fluctuations iCoh | ' num2str(num_boots) ' Runs'])
xlabel('Mean Estimates')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')

figure;
ksdensity(baseline_var_boot_icoh, 'Function', 'pdf');
hold on;
ksdensity(low_dose_var_boot_icoh, 'Function', 'pdf');
hold on;
ksdensity(high_dose_var_boot_icoh, 'Function', 'pdf');
title(['Boostrapped estimates of the Variance of Fluctuations iCoh | ' num2str(num_boots) ' Runs'])
xlabel('Variance Estimates')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')

%% 95th Percent Confidence Intervals on Boostrapped Estimates of the Variance of Fluctuations

ConfInt = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % calculate confidence interval on bootstrapped estimates 
pct = 95; % set percent CI

baseline_var_boot_psi_CI = ConfInt(baseline_var_boot_psi, pct);
low_dose_var_boot_psi_CI = ConfInt(low_dose_var_boot_psi, pct);
high_dose_var_boot_psi_CI = ConfInt(high_dose_var_boot_psi, pct);

baseline_var_boot_icoh_CI = ConfInt(baseline_var_boot_icoh, pct);
low_dose_var_boot_icoh_CI = ConfInt(low_dose_var_boot_icoh, pct);
high_dose_var_boot_icoh_CI = ConfInt(high_dose_var_boot_icoh, pct);


%% Variance PDFs with reported stats

a = 0.4; % set alpha, or, transparency for line fill

figure;
[f, x] = ksdensity(baseline_var_boot_psi, 'Function', 'pdf');
fill(x, f, [0 0.4470 0.7410], 'facealpha', a, 'LineStyle','none')
hold on;
[f, x] = ksdensity(low_dose_var_boot_psi, 'Function', 'pdf');
fill(x, f, [0.8500 0.3250 0.0980], 'facealpha', a, 'LineStyle','none')
hold on;
[f, x] = ksdensity(high_dose_var_boot_psi, 'Function', 'pdf');
fill(x, f, [0.4940 0.1840 0.5560], 'facealpha', a, 'LineStyle','none')
title(['Boostrapped Variance of Fluctuations PSI | ' num2str(num_boots) ' Runs'])
xlabel('Variance Estimates')
ylabel('Density')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')
xline(baseline_var_boot_psi_CI(1), 'b','HandleVisibility','off')
xline(baseline_var_boot_psi_CI(2), 'b','HandleVisibility','off')
xline(low_dose_var_boot_psi_CI(1), 'r','HandleVisibility','off')
xline(low_dose_var_boot_psi_CI(2), 'r','HandleVisibility','off')
xline(high_dose_var_boot_psi_CI(1), 'm','HandleVisibility','off')
xline(high_dose_var_boot_psi_CI(2), 'm','HandleVisibility','off')


figure;
[f, x] = ksdensity(baseline_var_boot_icoh, 'Function', 'pdf');
fill(x, f, [0 0.4470 0.7410], 'facealpha', a, 'LineStyle','none')
hold on;
[f, x] = ksdensity(low_dose_var_boot_icoh, 'Function', 'pdf');
fill(x, f, [0.8500 0.3250 0.0980], 'facealpha', a, 'LineStyle','none')
hold on;
[f, x] = ksdensity(high_dose_var_boot_icoh, 'Function', 'pdf');
fill(x, f, [0.4940 0.1840 0.5560], 'facealpha', a, 'LineStyle','none')
title(['Boostrapped Variance of Fluctuations iCoh | ' num2str(num_boots) ' Runs'])
xlabel('Variance Estimates')
ylabel('Density')
legend('Baseline', '0.2 \mug/mL Ketamine', '0.4 \mug/mL Ketamine')
xline(baseline_var_boot_icoh_CI(1), 'b','HandleVisibility','off')
xline(baseline_var_boot_icoh_CI(2), 'b','HandleVisibility','off')
xline(low_dose_var_boot_icoh_CI(1), 'r','HandleVisibility','off')
xline(low_dose_var_boot_icoh_CI(2), 'r','HandleVisibility','off')
xline(high_dose_var_boot_icoh_CI(1), 'm','HandleVisibility','off')
xline(high_dose_var_boot_icoh_CI(2), 'm','HandleVisibility','off')

%% Export Data To CSV for Python Plotting

all_flux = [BaselineFluctuations_psi(:), lowKFluctuations_psi(:), highKFluctuations_psi(:), BaselineFluctuations_icoh(:), lowKFluctuations_icoh(:), highKFluctuations_icoh(:)];
writematrix(all_flux,'C:\Users\diego\Desktop\all_fluctuations.csv'); 

all_boots = [baseline_var_boot_psi, low_dose_var_boot_psi, high_dose_var_boot_psi, baseline_var_boot_icoh, low_dose_var_boot_icoh, high_dose_var_boot_icoh];
writematrix(all_boots,'C:\Users\diego\Desktop\all_bootstrap_var.csv'); 

save('C:\Users\diego\Desktop\all_matrices.mat', 'Penn_baseline_alpha_psi', 'Penn_low_dose_alpha_psi', 'Penn_high_dose_alpha_psi', 'Penn_baseline_iCoh_segmented', 'Penn_low_dose_iCoh_segmented', 'Penn_high_dose_iCoh_segmented')

writematrix(highKFluctuations_psi, 'C:\Users\diego\Desktop\high_psi_fluctuations.csv')
writematrix(lowKFluctuations_psi, 'C:\Users\diego\Desktop\low_psi_fluctuations.csv')
writematrix(BaselineFluctuations_psi, 'C:\Users\diego\Desktop\baseline_psi_fluctuations.csv')


%% Pitman-Morgan Test for eqauality of variances

[h1,p1,ratio1] = PitmanMorganTest(BaselineFluctuations_psi(:), lowKFluctuations_psi(:));
[h2,p2,ratio2] = PitmanMorganTest(BaselineFluctuations_psi(:), highKFluctuations_psi(:));
[h3,p3,ratio3] = PitmanMorganTest(lowKFluctuations_psi(:), highKFluctuations_psi(:));

[h1c,p1c,ratio1c] = PitmanMorganTest(BaselineFluctuations_icoh(:), lowKFluctuations_icoh(:));
[h2c,p2c,ratio2c] = PitmanMorganTest(BaselineFluctuations_icoh(:), highKFluctuations_icoh(:));
[h3c,p3c,ratio3c] = PitmanMorganTest(lowKFluctuations_icoh(:), highKFluctuations_icoh(:));

T_PSI = table([h1; h2; h3], [p1; p2; p3], [ratio1; ratio2; ratio3], 'VariableNames',{'Is Stat. Sig.?','p-val', 'Var. Ratio'}, 'RowNames',{'Baseline-LowDose','Baseline-HighDose','LowDose-HighDose'});
T_iCOH = table([h1c; h2c; h3c], [p1c; p2c; p3c], [ratio1c; ratio2c; ratio3c], 'VariableNames',{'Is Stat. Sig.?','p-val', 'Var. Ratio'}, 'RowNames',{'Baseline-LowDose','Baseline-HighDose','LowDose-HighDose'});

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh([p1, p2, p3]);
[hx, crit_px, adj_ci_cvrgx, adj_px] = fdr_bh([p1c, p2c, p3c]);
