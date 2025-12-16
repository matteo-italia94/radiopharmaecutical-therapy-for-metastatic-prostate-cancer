% main_simulation.m
%
% Author: Matteo Italia, ORCID ID 0000-0001-9633-8396 (email: matteo.italia@uclm.es)
%
% Main script to run in silico radiopharmaceutical therapy (RPT) trials.
% - Runs N virtual patients through a treatment protocol
% - Calls simulate_virtual_patient_under_RPT function (which calls RPT_dynamics_virtual_patient function) for per-patient ODE simulation and biologically effective dose (BED), absorbed doses (AD), and overall survival (OS) outputs
% - Collects aggregated statistics and produces basic plots / survival analysis
%
% Dependencies: MATLAB (recommend R2019b+), Statistics / Parallel Toolbox (parfor).

%clear all; close all; clc;
tic % for recording simulation time
rng(1); % for reproducibility
%% ------------------------
%% User / experiment settings
%% ------------------------
recording_all = 0;    % store full time series per patient (memory heavy)
different_patient = 0;  % if 0 -> use the same virtual patients, i.e., the same parameters (params_VPs_vector must exist in the workspace or be loaded)
				    % if 1-> generate new virtual patients

if ~exist('parameters_500_VPs.mat','file')
    fprintf('Missing file: parameters_500_VPs.mat. Please add it to the repository to reproduce Results presented in the manuscript.');
end
load('parameters_500_VPs.mat'); % to reproduce Results presented in the article


%% Treatment settings (e.g., standard treatment: 6 injections, each injection of 7.4 GBq and separated by 6-week one from the consecutive)
n_inj = 6; % number of injections
tot_activity = 44.4;   % total administered activity across all cycles [GBq]
tot_activity = tot_activity*1e3;   % total administered activity across all cycles [MBq]
dose = tot_activity / n_inj; % injection dose  7.4 1e3 [MBq] 
cycle_dose_vector = repmat(dose, n_inj, 1);  % vector of per-cycle activities [MBq]
%cycle_dose_vector = [7.4;37]*1e3; % vector of per-cycle activities for 2* regimen[MBq]
cycle_leght = 6; % 6-weeks
cycle_length_vector = repmat(cycle_leght, n_inj, 1) * 7 * 24; % vector of cycle lengths: 6 weeks per cycle -> convert to hours [h] 


%% Simulation parameters
N_patients = 500;    % number of virtual patients (VPs) in the in silico trial (virtual cohort size)
N_param = 30;      % number of parameters per virtual patient (used in params_VPs_vector)
delta_t_sim = 0.1;   % integration time step [hours]

% Physical constants / common parameters
Ee_Lu177 = 0.0849;   % energy emitted per beta decay [mJ·MBq−1·h−1]
S_gamma_WB = 0.00185*1e-3; % Self-dose S-value for whole body (Gy·Kg^(2/3)·MBq−1·h−1)
lambda_phy = 0.0044;  % physical decay rate (h^-1) of 177Lu 

%% Initialize storage
if different_patient
  params_VPs_vector = zeros(N_patients, N_param); % clear vector and generate new virtual patients
end

% vectors to collect the in silico trial results
T_survival_day = zeros(N_patients,1);  % OS [day]
alive = zeros(N_patients,1); % vector indicating alive VPs at the end of the trial (->1 for alive VPs)
death = zeros(N_patients,1); % vector indicating dead VPs at the end of the trial (->1 for dead VPs)
Event = cell(N_patients,1);
Times_KM = zeros(N_patients,1);
outcomes = zeros(N_patients,5); % [Times_KM, BED_kidney, BED_PG, BED_LG, BED_RM]

% Per-patient aggregated outputs (preallocate)
AD_all_PG = zeros(N_patients,1); % ADs in the parotid glands
AD_all_LG = zeros(N_patients,1); % ADs in the lacrimal glands
AD_all_RM = zeros(N_patients,1); % AD in the red marrow
AD_all_K = zeros(N_patients,1); % AD in the kidneys
AD_all_T = zeros(N_patients,1); % AD in tumor
AD_all_WB = zeros(N_patients,1); % AD in the whole body


BED_Kidneys = zeros(N_patients,1); % BED in the kidneys
BED_lacrimal = zeros(N_patients,1); % BED in the lacrimal glands
BED_parotid = zeros(N_patients,1); % BED in the parotid glands
BED_redmarrow = zeros(N_patients,1); % BED in the red marrow

%% ------------------------
%% Generate virtual patient parameter vectors
% different_patient must be 1 to generate virtual patients 
%% ------------------------
if different_patient
  for i = 1:N_patients
    % Draw physiological and radiobiological parameters from uniform/random distributions. 
    M_TP0 = 100 * rand; % initial proliferative tumor mass [g] 
    params_VPs_vector(i,1) = M_TP0;

    MK = (231 + (503 - 231) * rand); % total kidneys mass [g] 
    params_VPs_vector(i,2) = MK;

    MPG = 2*(15 + rand*(30-15)); % parotid glands mass [g]
    params_VPs_vector(i,3) = MPG;

    MLG = 2*(0.5 + (0.780 - 0.5) * rand) ; % lacrimal glands mass [g]
    params_VPs_vector(i,4) = MLG;

    MRM = (1.372 + (2.852 - 1.372) * rand) * 1.03 * 1e3; % red marrow mass [g] using red marrow density = 1.03 g/cm³

    params_VPs_vector(i,5) = MRM;

    alpha = 0.01 + (0.15 - 0.01) * rand; % linear coefficient in the linear-quadratic model [Gy-1] 
    params_VPs_vector(i,6) = alpha;

    betha = 0.003 + (0.058 - 0.003) * rand; % quadratic coefficient in the linear-quadratic model [Gy-2] 
    params_VPs_vector(i,7) = betha;

    epsilon = rand(); % sublethal lesions interaction probability [non-dimensional]
    p = sqrt(betha / epsilon);
    params_VPs_vector(i,8) = p; % yield of sublethal lesions per unit AD [Gy-1]
    params_VPs_vector(i,9) = epsilon;

    muT = 0.3648; % sublethal damage repair rate in tumor [h-1]
    params_VPs_vector(i,10) = muT;

    kc = (1 + (10 - 1) * rand) * 1e-4; % tumor clearance rate [h-1]
    params_VPs_vector(i,11) = kc;

    kg = (1 + (19 - 1) * rand) * 1e-4; % tumor growth rate [h-1]
    params_VPs_vector(i,12) = kg;

    params_VPs_vector(i,13) = lambda_phy; % physical decay rate of 177Lu [h^-1]

    lambda_r_T = 1e-5 + (0.045 - 1e-5) * rand; % Biological Clearance Rate for the tumor [h^-1] 
    params_VPs_vector(i,14) = lambda_r_T;

    lambda_r_K = 0.004 + (0.032 - 0.004) * rand; % Biological Clearance Rate for the kidneys [h^-1]  
    params_VPs_vector(i,15) = lambda_r_K;

    lambda_r_PG = 0.0117 + (0.0303 - 0.0117) * rand; % Biological Clearance Rate for the parotid glands [h^-1] 
    params_VPs_vector(i,16) = lambda_r_PG;

    lambda_r_LG = 0.0204; % % Biological Clearance Rate for the lacrimal glands [h^-1] 
    params_VPs_vector(i,17) = lambda_r_LG;

    lambda_r_RM = 0.006 + (0.0438 - 0.006) * rand; % Biological Clearance Rate for the red marrow [h^-1] 
    params_VPs_vector(i,18) = lambda_r_RM;

    lambda_r_RoB = 0.0032 + (0.0286 - 0.0032) * rand; % Biological Clearance Rate for the whole body [h^-1] 
    params_VPs_vector(i,19) = lambda_r_RoB;

    alpha_over_beta_K = 1.4 + (4.3 - 1.4) * rand; % Biological Clearance Rate for the tumor [h^-1] 
    params_VPs_vector(i,20) = alpha_over_beta_K; 

    mu_K = 0.1 + (0.57 - 0.1) * rand; % Radiobiological α/β for kidneys [Gy]
    params_VPs_vector(i,21) = mu_K;

    MWB = (60 + (100-60) * rand)*1e3; % whole-body mass [g]
    params_VPs_vector(i,22) = MWB;

    M_TDeath = 1000; % death tumor mass threshold [g]
    params_VPs_vector(i,23) = M_TDeath;

    % Standardize Uptake Values (SUV)
	SUV_mean_T = 3 + 77*rand;     % tumor SUV mean [non-dimensional]   
    params_VPs_vector(i,24) = SUV_mean_T;

    SUV_mean_PG = 5 + (15-5) * rand; % parotid glands SUV mean [non-dimensional]
    SUV_mean_LG = 1.5 + (7-1.5) * rand; % lacrimal glands SUV mean [non-dimensional]
    SUV_mean_K = 4 + (30-4) * rand; % kidneys SUV mean [non-dimensional]
    SUV_mean_BM = 0.3 + (1.4-0.3) * rand; % bone marrow SUV mean [non-dimensional]
    MRoB = MWB-MRM-MK-MPG-MLG-M_TP0; % mass of the rest of the body [g]
    SUV_mean_RoB=-1; % rest of the body SUV mean [non-dimensional]
    while SUV_mean_RoB<0 % check biology consistency (SUV mean must be non negative)
            SUV_mean_WB = 1; % whole body SUV mean [non-dimensional]
            SUV_mean_RoB = 1/MRoB * (SUV_mean_WB*MWB - SUV_mean_BM*MRM - SUV_mean_PG*MPG - SUV_mean_LG*MLG - SUV_mean_K*MK - SUV_mean_T*M_TP0);
    end    

    params_VPs_vector(i,25) = SUV_mean_BM;
    params_VPs_vector(i,26) = SUV_mean_K;
    params_VPs_vector(i,27) = SUV_mean_LG;
    params_VPs_vector(i,28) = SUV_mean_PG;
    params_VPs_vector(i,29) = SUV_mean_RoB;  
    params_VPs_vector(i,30) = MRoB;  

  end
end

%% ------------------------
%% Time vector for full simulation horizon
%% ------------------------
t_50months = 50 * 4 * 7 * 24;  % 50 months in hours [h] 
t_sim = 0:delta_t_sim:t_50months; % simulation time vector [hours]

%% If recording_all, preallocate large matrices (memory heavy)
if recording_all
  Tumor_Mass_All = zeros(N_patients, length(t_sim));
  Activity_K_All = zeros(N_patients, length(t_sim));
  Activity_RM_All = zeros(N_patients, length(t_sim));
  Activity_LG_All = zeros(N_patients, length(t_sim));
  Activity_PG_All = zeros(N_patients, length(t_sim));
  Activity_WB_All = zeros(N_patients, length(t_sim));
  Activity_T_All = zeros(N_patients, length(t_sim));
  AD_T_All = zeros(N_patients, length(t_sim));
end

%% ------------------------
%% Run per-patient simulations (parallel loop recommended)
%% ------------------------
% NOTE: the function simulate_virtual_patient_under_RPT implement treatment applications, 
% using the function RPT_dynamics_virtual_patient to implement the ODE integration, and 
% return the patient survival, the time series of activities, tumor compartments, BEDs, etc.
parfor pat = 1:N_patients

[T_survival_day(pat), death(pat), alive(pat), ...
  Event{pat}, Times_KM(pat), ...
  M_TP, M_TSD, M_TD, A_kidney, A_marrow, A_lacrimal, A_parotid, ...
  A_rest_body, A_tumor, AD_tumor, ...
  BED_K, BED_RM, BED_LG, BED_PG] = ...
    simulate_virtual_patient_under_RPT(pat, params_VPs_vector, cycle_dose_vector, cycle_length_vector, delta_t_sim, n_inj, Ee_Lu177, S_gamma_WB);

  M_TT = M_TP + M_TSD + M_TD; % The total tumor mass (M_TT) [g]. It is the sum of proliferating tumor cells (M_TP), sublethally damaged cells (M_TSD), and lethally damaged cells (M_TD)

    % store single VP detailed results to calculate virtual cohort statistics
  if recording_all 
    Tumor_Mass_All(pat, :) = M_TT;
    Activity_K_All(pat, :) = A_kidney;
    Activity_RM_All(pat, :) = A_marrow;
    Activity_LG_All(pat, :) = A_lacrimal;
    Activity_PG_All(pat, :) = A_parotid;
    Activity_WB_All(pat, :) = A_rest_body;
    Activity_T_All(pat, :) = A_tumor;
    AD_T_All(pat, :) = AD_tumor;
  end

  % Aggregate end-of-simulation metrics (trapezoidal integration)
  MK = params_VPs_vector(pat,2);
  MPG = params_VPs_vector(pat,3);
  MLG = params_VPs_vector(pat,4);
  MRM = params_VPs_vector(pat,5);
  MWB = params_VPs_vector(pat,22);
  MRoB = params_VPs_vector(pat,30);  

  % Parotid dose
  SPG = Ee_Lu177 / MPG;
  A_cum_PG = trapz(t_sim, A_parotid);
  DPG = SPG * A_cum_PG;

  % Lacrimal dose
  SLG = Ee_Lu177 / MLG;
  A_cum_LG = trapz(t_sim, A_lacrimal);
  DLG = SLG * A_cum_LG;

  % Kidney dose
  SK = Ee_Lu177 / MK;
  A_cum_K = trapz(t_sim, A_kidney);
  DK = SK * A_cum_K;

  % Rest-of-the-body dose
  SRoB = Ee_Lu177 / MRoB;
  A_cum_RoB = trapz(t_sim, A_rest_body);
  DWB = SRoB * A_cum_RoB;

  % Red marrow dose: plus gamma contribution from whole-body 
  SRM = Ee_Lu177 / MRM;
  A_cum_RM = trapz(t_sim, A_marrow);
  A_cum_T = trapz(t_sim, A_tumor);
  A_cum_WB=A_cum_RM+A_cum_PG+A_cum_LG+A_cum_K+A_cum_T+A_cum_RoB;
  DRM = SRM * A_cum_RM + S_gamma_WB/((MWB*1e-3)^(2/3)) * A_cum_WB;

  % Save aggregated metrics
  AD_all_PG(pat) = DPG;
  AD_all_LG(pat) = DLG;
  AD_all_RM(pat) = DRM;
  AD_all_K(pat) = DK;
  AD_all_WB(pat) = DWB;
  AD_all_T(pat) = AD_tumor(end);
  BED_Kidneys(pat) = BED_K;
  BED_parotid(pat) = BED_PG;
  BED_lacrimal(pat) = BED_LG;
  BED_redmarrow(pat) = BED_RM;

  outcomes(pat,:) = [Times_KM(pat), BED_K, BED_PG, BED_LG, BED_RM];
end

%% ------------------------
%% Toxicity thresholds   
%% ------------------------
max_BED_K = 39;  % Kidney BED threshold [Gy] 
max_BED_LG = 68; % Lacrimal BED threshold [Gy] 
max_BED_PG = 31.5;  % Parotid BED threshold [Gy] 
max_BED_RM = 2.02;  % Red marrow BED threshold [Gy] 

%% Uncomment following lines to compute max BED for parotid glands (PG) and red marrow (RM) using original formulas
% alpha_over_beta_PG = 0.053 / 0.0118;
% mu_PG = 0.0077 * 60; % convert min^-1 to h^-1 if original was per minute
% lambda_r_PG_med = 0.0117 + (0.0303 - 0.0117) * 0.5;
% T_rep_PG = log(2) / mu_PG;
% lambda_eff_PG = lambda_r_PG_med + lambda_phy;
% T_eff_PG = log(2) / lambda_eff_PG;
% max_D_PG=26; %[Gy]
% max_BED_PG = max_D_PG + 1/alpha_over_beta_PG * (T_rep_PG/(T_rep_PG * T_eff_PG)) * max_D_PG^2;

%alpha_over_beta_RM = 15;
%mu_RM = log(2) / 1.5;
%T_rep_RM = log(2) / mu_RM;
%lambda_r_RM_med = 0.0074 + (0.077 - 0.0074) * 0.5;
%lambda_eff_RM = lambda_r_RM_med + lambda_phy;
%T_eff_RM = log(2) / lambda_eff_RM;
%max_D_RM = 2; %[Gy]
%max_BED_RM = max_D_RM + 1/alpha_over_beta_RM * (T_rep_RM/(T_rep_RM * T_eff_RM)) * max_D_RM^2;

%% ------------------------
%% Compute statistics and toxicity probabilities
%% ------------------------
statistics = 1;
if statistics
  fprintf('Number of injections: %d\n', n_inj);
  fprintf('Cycle length (week): %.3f\n', cycle_leght);
  fprintf('Dose per injection (MBq): %.3f\n', dose);

  BED_PG_all = [mean(BED_parotid), std(BED_parotid)];
  prob_toxicity_BED_PG = sum(BED_parotid > max_BED_PG) / N_patients;

  BED_LG_all = [mean(BED_lacrimal), std(BED_lacrimal)];
  prob_toxicity_BED_LG = sum(BED_lacrimal > max_BED_LG) / N_patients;

  BED_RM_all = [mean(BED_redmarrow), std(BED_redmarrow)];
  prob_toxicity_BED_RM = sum(BED_redmarrow > max_BED_RM) / N_patients;

  BED_K_all = [mean(BED_Kidneys), std(BED_Kidneys)];
  prob_toxicity_BED_K = sum(BED_Kidneys > max_BED_K) / N_patients;

  AD_WB = [mean(AD_all_WB), std(AD_all_WB)];
  AD_T = [mean(AD_all_T), std(AD_all_T)];

  % Print a compact summary
  fprintf('Prob toxicity (kidney, BED): %.3f\n', prob_toxicity_BED_K);
  fprintf('Prob toxicity (red marrow, BED): %.3f\n', prob_toxicity_BED_RM);
  fprintf('Prob toxicity (parotid, BED): %.3f\n', prob_toxicity_BED_PG);
  fprintf('Prob toxicity (lacrimal, BED): %.3f\n', prob_toxicity_BED_LG);
end


%% ------------------------
%% Optional: plotting when recording_all is true
%% ------------------------
if recording_all
  % compute percentiles / median and show tumor mass dynamics
  median_Tumor_Mass = median(Tumor_Mass_All, 1);
  CI_90_low = prctile(Tumor_Mass_All, 5, 1);
  CI_50_low = prctile(Tumor_Mass_All, 25, 1);
  CI_50_high = prctile(Tumor_Mass_All, 75, 1);
  CI_90_high = prctile(Tumor_Mass_All, 95, 1);
  t_weeks = t_sim / (24 * 7);
  figure; hold on;
  plot(t_weeks, median_Tumor_Mass, 'k', 'LineWidth', 2);
  plot(t_weeks, CI_50_high, 'b--', 'LineWidth', 1.5);
  plot(t_weeks, CI_50_low, 'b--', 'LineWidth', 1.5);
  plot(t_weeks, CI_90_high, 'r--', 'LineWidth', 1);
  plot(t_weeks, CI_90_low, 'r--', 'LineWidth', 1);
  xlabel('Time (weeks)');
  ylabel('Total Tumor Mass (g)');
  title(sprintf('Tumor Growth Dynamics in Virtual In Silico Trial (N=%d)', N_patients));
  legend({'Median', '75% quartile', '25% quartile','95% quartile', '5% quartile'}, 'Location', 'northeast');
  grid on; xlim([0 max(t_weeks)]);
end

%% ------------------------
%% Survival analysis 
%% ------------------------
EventNum = strcmp(Event, 'Dead');
[f, x] = ecdf(Times_KM, 'Censoring', EventNum == 0);
S = 1 - f;
% Find the median (time at which survival drops below 0.5)
index = find(S <= 0.5, 1, 'first');

if isempty(index)
    median_OS = NaN; % median not reached
    disp('Median OS not reached (less than 50% events)');
else
    median_OS = x(index);
    fprintf('Median OS = %.3f (months)\n', median_OS);
end


toc