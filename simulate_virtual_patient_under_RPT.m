function [T_survival_day, death, alive, Event, Times_KM, M_TP, M_TSD, M_TD, A_kidney, A_marrow, A_lacrimal,A_parotid, A_rest_body, A_tumor, AD_tumor, BED_K, BED_RM, BED_LG, BED_PG] ...
         = simulate_virtual_patient_under_RPT( ...
           pat, params_VPs_vector, cycle_dose_vector, ...
           cycle_length_vector, delta_t_sim, NC, Ee_Lu177, S_gamma_WB)
%--------------------------------------------------------------------------
% simulate_virtual_patient_under_RPT
%
% Author: Matteo Italia, ORCID ID 0000-0001-9633-8396 (email: matteo.italia@uclm.es)
%
% Simulates a virtual patient (VP)’s tumor and organ kinetics under
% radiopharmaceutical therapy (RPT) with 177Lu-labeled PSMA.
%
% The model integrates tumor growth dynamics, multi-organ activity kinetics,
% absorbed dose (AD) and biologically effective dose (BED) accumulation,
% and survival time estimation.
%
% The function calls the RPT_dynamics_virtual_patient function to simulate tumor and activities evolutions between consecutive injections.
%
% INPUTS:
%   pat     - patient index
%   params_VPs_vector   - matrix of virtual patient parameters
%   cycle_dose_vector   - administered activity per cycle [MBq]
%   cycle_length_vector - length of each treatment cycle [h]
%   delta_t_sim   - simulation time step [h]
%   recording_all    - boolean, record all activities and tumor mass
%   NC      - number of treatment cycles
%   Ee_Lu177   - mean beta energy emitted per decay [Gy·kg/GBq·h]
%   S_gamma_WB - Self-dose S-value for whole body [Gy·Kg^(2/3)·MBq−1·h−1]
%
% OUTPUTS:
%   T_survival_day   - simulated survival time [days]
%   death/alive   - binary indicators of death/survival
%   Event      - 'Dead' or 'Living'
%   Times_KM   - time to event in months (for Kaplan-Meier)
%   BEDK    - kidney biologically effective dose (BED)
%   BED_RM, BED_LG, BED_PG, BED_K - BED for red marrow, lacrimal, parotid glands, and kidneys, respectively
%
%--------------------------------------------------------------------------

%% Constants and organ radiobiological parameters
% Parotid gland (PG)
alpha_over_beta_PG = 0.053 / 0.0118; % (Gy)
mu_PG = 0.0077 * 60; % (1/h)

% Red marrow (RM)
alpha_over_beta_RM = 15; % (Gy)
mu_RM = log(2) / 1.5; % (1/h)

% Lacrimal gland (LG)
T_half_rep_LG = 1.5; % (h)
alpha_over_beta_LG = 1 / 4.5;% (Gy)

%% Extract patient-specific parameters
M_TP0 = params_VPs_vector(pat, 1);
MK = params_VPs_vector(pat, 2);
MPG = params_VPs_vector(pat, 3);
MLG = params_VPs_vector(pat, 4);
MRM = params_VPs_vector(pat, 5);

lambda_phy= params_VPs_vector(pat, 13);
lambda_r_T= params_VPs_vector(pat, 14);
lambda_r_K = params_VPs_vector(pat, 15);
lambda_r_PG = params_VPs_vector(pat, 16);
lambda_r_LG = params_VPs_vector(pat, 17);
lambda_r_RM = params_VPs_vector(pat, 18);
lambda_r_WB = params_VPs_vector(pat, 19);

alpha_over_beta_K = params_VPs_vector(pat, 20);
mu_K = params_VPs_vector(pat, 21);

MWB = params_VPs_vector(pat, 22);
death_size = params_VPs_vector(pat, 23);
SUV_mean_T = params_VPs_vector(pat, 24);
SUV_mean_RM = params_VPs_vector(pat, 25);
SUV_mean_K = params_VPs_vector(pat, 26);
SUV_mean_LG = params_VPs_vector(pat, 27);
SUV_mean_PG = params_VPs_vector(pat, 28);
SUV_mean_RoB = params_VPs_vector(pat, 29);
MRoB = params_VPs_vector(pat, 30);


%% Initialize tumor and organ compartments
M_TSD0 = 0; 
M_TD0 = 0;
M_TT = M_TP0 + M_TD0 + M_TSD0;
denominator=SUV_mean_T * M_TT+SUV_mean_PG * MPG+SUV_mean_LG * MLG+SUV_mean_K * MK+SUV_mean_RM * MRM+SUV_mean_RoB * MRoB;
fT = SUV_mean_T * M_TT / denominator;
fK = SUV_mean_K * MK / denominator;
fLG = SUV_mean_LG * MLG / denominator;
fPG = SUV_mean_PG * MPG / denominator;
fRM = SUV_mean_RM * MRM / denominator;
fRoB = SUV_mean_RoB * MRoB / denominator;
AWB = cycle_dose_vector(1);
A_T = fT * AWB;
A_kidney0 = fK * AWB;
A_parotid0 = fPG * AWB;
A_lacrimal0 = fLG * AWB;
A_marrow0 = fRM * AWB;
A_rest_body0 = fRoB*AWB;
AD_tumor0 = 0;

%% Initial conditions vector
y = [M_TP0, M_TSD0, M_TD0, A_T, A_kidney0, A_parotid0, A_lacrimal0, A_rest_body0, A_marrow0, AD_tumor0];

num_vars = numel(y);
options = odeset('NonNegative', 1:num_vars);
options = odeset(options, 'RelTol',1e-8,'AbsTol',1e-10);
params_pat = params_VPs_vector(pat, :);

%% --- Simulate first treatment cycle ---
tcycle = 0:delta_t_sim:cycle_length_vector(1);
[t, Y] = ode45(@(t,y) RPT_dynamics_virtual_patient(t, y, params_pat, Ee_Lu177), tcycle, y, options);

[M_TP, M_TSD, M_TD, A_tumor, A_kidney, A_parotid, A_lacrimal, A_rest_body, A_marrow, AD_tumor] = deal(Y(:,1), Y(:,2), Y(:,3), Y(:,4), Y(:,5), Y(:,6), Y(:,7), Y(:,8), Y(:,9), Y(:,10));

Time = t;

%% --- Initialize dose and BED accumulators ---
BED_K = 0; BED_RM = 0; BED_LG = 0; BED_PG = 0;
% AD1_new_T = 0; AD2_new_T = 0;

% Compute absorbed and biologically effective dose for first cycle

%kidney
SK = Ee_Lu177 / MK;
A_cum_K = trapz(Time, A_kidney);
DK_i = SK * A_cum_K; %dose to kidneys
lambda_eff_K_pat = lambda_r_K + lambda_phy; % (h^-1)
%NOTE: mu_K = log(2)/T_0_5_rep;
% and T 1/2_eff = ln(2)/ λ eff;
T_half_rep_K=log(2)/mu_K;
T_half_eff_K=log(2)/lambda_eff_K_pat;
BEDK_i=DK_i+(1/alpha_over_beta_K)*T_half_rep_K/(T_half_rep_K+T_half_eff_K)*DK_i^2;
BED_K=BED_K+BEDK_i;

%Parotid gland
SPG = Ee_Lu177 / MPG;
A_cum_PG = trapz(Time, A_parotid);
DPG_i = SPG * A_cum_PG; %dose to pg
lambda_eff_PG_pat= lambda_r_PG + lambda_phy; % (h^-1)
T_half_rep_PG=log(2)/mu_PG;
T_half_eff_PG=log(2)/lambda_eff_PG_pat;
BED_PG_i=DPG_i+(1/alpha_over_beta_PG)*T_half_rep_PG/(T_half_rep_PG+T_half_eff_PG)*DPG_i^2;
BED_PG=BED_PG+BED_PG_i;

%lacrimal gland
SLG = Ee_Lu177 / MLG;
A_cum_LG = trapz(Time, A_lacrimal);
DLG_i = SLG * A_cum_LG; %dose to lg
lambda_eff_LG_pat= lambda_r_LG + lambda_phy; % (h^-1)
T_half_eff_LG=log(2)/lambda_eff_LG_pat;
BED_LG_i=DLG_i+(1/alpha_over_beta_LG)*T_half_rep_LG/(T_half_rep_LG+T_half_eff_LG)*DLG_i^2;
BED_LG=BED_LG+BED_LG_i;

%red marrow
SRM = Ee_Lu177 / MRM;
A_cum_RM = trapz(Time, A_marrow);
A_cum_T = trapz(Time, A_tumor);
A_cum_RoB = trapz(Time, A_rest_body);
A_cum_WB=A_cum_RM+A_cum_PG+A_cum_LG+A_cum_K+A_cum_RoB+A_cum_T; %sum all A_cum
DRM_i = SRM * A_cum_RM+ S_gamma_WB/((MWB*1e-3)^(2/3))*A_cum_WB; %dose to RM
lambda_eff_RM_pat= lambda_r_RM + lambda_phy; % (h^-1)
T_half_rep_RM=log(2)/mu_RM;
T_half_eff_RM=log(2)/lambda_eff_RM_pat;
BED_RM_i=DRM_i+(1/alpha_over_beta_RM)*T_half_rep_RM/(T_half_rep_RM+T_half_eff_RM)*DRM_i^2;
BED_RM=BED_RM+BED_RM_i;

%% --- Simulate subsequent cycles ---
for cycle = 2:NC
    % Update tumor activity fraction
    M_TT = M_TP(end) + M_TSD(end) + M_TD(end);    
    % Update activities for next cycle
    AWB = cycle_dose_vector(cycle);
    
    % update fractions 
    denominator=SUV_mean_T * M_TT+SUV_mean_PG * MPG+SUV_mean_LG * MLG+SUV_mean_K * MK+SUV_mean_RM * MRM+SUV_mean_RoB * MRoB;
    fT = SUV_mean_T * M_TT/ denominator;
    fK = SUV_mean_K * MK / denominator;
    fLG = SUV_mean_LG * MLG / denominator;
    fPG = SUV_mean_PG * MPG / denominator;
    fRM = SUV_mean_RM * MRM / denominator;
    fRoB = SUV_mean_RoB * MRoB / denominator;
    %calculate activity for each compartment
    A_T = fT * AWB;
    AK = fK * AWB;
    APG = fPG * AWB;
    ALG = fLG * AWB;
    ARM = fRM * AWB;
    ARoB = fRoB*AWB;
    % Initial conditions for next cycle
    y = [M_TP(end), M_TSD(end), M_TD(end), A_tumor(end) + A_T, A_kidney(end) + AK, A_parotid(end) + APG, A_lacrimal(end) + ALG, A_rest_body(end) + ARoB, A_marrow(end) + ARM, AD_tumor(end)];
    t_start = sum(cycle_length_vector(1:cycle-1));
    t_end = t_start + cycle_length_vector(cycle);
    tspan = t_start:delta_t_sim:t_end;
    [t, Y] = ode45(@(t,y) RPT_dynamics_virtual_patient(t, y, params_pat, Ee_Lu177),tspan, y, options);

     % Append results
     M_TP = [M_TP; Y(2:end,1)];
     M_TSD = [M_TSD; Y(2:end,2)];
     M_TD = [M_TD; Y(2:end,3)];
     A_tumor = [A_tumor; Y(2:end,4)];
     A_kidney = [A_kidney; Y(2:end,5)];
     A_parotid = [A_parotid; Y(2:end,6)];
     A_lacrimal = [A_lacrimal; Y(2:end,7)];
     A_rest_body = [A_rest_body; Y(2:end,8)];
     A_marrow = [A_marrow; Y(2:end,9)];
     AD_tumor = [AD_tumor; Y(2:end,10)];
     Time = [Time; t(2:end)];

    % Compute absorbed and BED for this cycle  
    %NOTE: use only the current cycle absorbed dose
    A_cum_K = trapz(t(2:end), Y(2:end,5));
    DK_i = SK * A_cum_K; %dose to kidneys        
    BEDK_i=DK_i+(1/alpha_over_beta_K)*T_half_rep_K/(T_half_rep_K+T_half_eff_K)*DK_i^2;
    BED_K=BED_K+BEDK_i;
    
    %Parotid gland
    A_cum_PG = trapz(t(2:end), Y(2:end,6));
    DPG_i = SPG * A_cum_PG; %dose to kidneys
    BED_PG_i=DPG_i+(1/alpha_over_beta_PG)*T_half_rep_PG/(T_half_rep_PG+T_half_eff_PG)*DPG_i^2;
    BED_PG=BED_PG+BED_PG_i;
    
    %lacrimal gland
    A_cum_LG = trapz(t(2:end), Y(2:end,7));
    DLG_i = SLG * A_cum_LG; %dose to kidneys
    BED_LG_i=DLG_i+(1/alpha_over_beta_LG)*T_half_rep_LG/(T_half_rep_LG+T_half_eff_LG)*DLG_i^2;
    BED_LG=BED_LG+BED_LG_i;
    
    %red marrow
    A_cum_RM = trapz(t(2:end), Y(2:end,9));
    A_cum_RoB = trapz(t(2:end), Y(2:end,8));
    A_cum_T = trapz(t(2:end), Y(2:end,4));
    A_cum_WB=A_cum_RM+A_cum_PG+A_cum_LG+A_cum_K+A_cum_RoB+A_cum_T; %sum all A_cum
    DRM_i = SRM * A_cum_RM+ S_gamma_WB/((MWB*1e-3)^(2/3))*A_cum_WB; %dose to RM
    BED_RM_i=DRM_i+(1/alpha_over_beta_RM)*T_half_rep_RM/(T_half_rep_RM+T_half_eff_RM)*DRM_i^2;
    BED_RM=BED_RM+BED_RM_i;
end

%% --- Post-therapy simulation ---
y = [M_TP(end), M_TSD(end), M_TD(end), A_tumor(end), A_kidney(end), A_parotid(end), A_lacrimal(end), A_rest_body(end), A_marrow(end), AD_tumor(end)];

t_50months = 50 * 4 * 7 * 24;
tspan = Time(end):delta_t_sim:t_50months;

[t, Y] = ode45(@(t,y) RPT_dynamics_virtual_patient(t, y, params_pat, Ee_Lu177), tspan, y, options);

% Append final results
M_TP = [M_TP; Y(2:end,1)];
M_TSD = [M_TSD; Y(2:end,2)];
M_TD = [M_TD; Y(2:end,3)];
A_tumor = [A_tumor; Y(2:end,4)];
A_kidney = [A_kidney; Y(2:end,5)];
A_parotid = [A_parotid; Y(2:end,6)];
A_lacrimal = [A_lacrimal; Y(2:end,7)];
A_rest_body = [A_rest_body; Y(2:end,8)];
A_marrow = [A_marrow; Y(2:end,9)];
AD_tumor = [AD_tumor; Y(2:end,10)];
Time = [Time; t(2:end)];

    %% --- BED estimations for remaining activity
    %NOTE: use only the new absorbed dose
    % kidneys
    A_cum_K = trapz(t(2:end), Y(2:end,5));
    DK_i = SK * A_cum_K; %dose to kidneys       
    BEDK_i=DK_i+(1/alpha_over_beta_K)*T_half_rep_K/(T_half_rep_K+T_half_eff_K)*DK_i^2;
    BED_K=BED_K+BEDK_i;
    
    %Parotid gland
    A_cum_PG = trapz(t(2:end), Y(2:end,6));
    DPG_i = SPG * A_cum_PG; %dose to kidneys
    BED_PG_i=DPG_i+(1/alpha_over_beta_PG)*T_half_rep_PG/(T_half_rep_PG+T_half_eff_PG)*DPG_i^2;
    BED_PG=BED_PG+BED_PG_i;

    %lacrimal gland
    A_cum_LG = trapz(t(2:end), Y(2:end,7));
    DLG_i = SLG * A_cum_LG; %dose to kidneys
    BED_LG_i=DLG_i+(1/alpha_over_beta_LG)*T_half_rep_LG/(T_half_rep_LG+T_half_eff_LG)*DLG_i^2;
    BED_LG=BED_LG+BED_LG_i;

    %red marrow
    A_cum_RM = trapz(t(2:end), Y(2:end,9));
    A_cum_RoB = trapz(t(2:end), Y(2:end,8));
    A_cum_T = trapz(t(2:end), Y(2:end,4));
    A_cum_WB=A_cum_RM+A_cum_PG+A_cum_LG+A_cum_K+A_cum_RoB+A_cum_T; %sum all A_cum
    DRM_i = SRM * A_cum_RM+ S_gamma_WB/((MWB*1e-3)^(2/3))*A_cum_WB; %dose to RM
    BED_RM_i=DRM_i+(1/alpha_over_beta_RM)*T_half_rep_RM/(T_half_rep_RM+T_half_eff_RM)*DRM_i^2;
    BED_RM=BED_RM+BED_RM_i;

    %% --- Survival evaluation ---
    if max(M_TP + M_TSD + M_TD) > death_size % check mortaly condition
     %T_survival_day = find(M_TP + M_TSD + M_TD > death_size, 1) / ((1 / delta_t_sim) * 24);
      idx = find(M_TP + M_TSD + M_TD > death_size,1);
      T_survival_day = Time(idx) / 24;
     death = 1; alive = 0; Event = 'Dead';
     Times_KM = T_survival_day / 28; % [months]
    else
     T_survival_day = 50 * 4 * 7;
     alive = 1; death = 0; Event = 'Living';
     Times_KM = 50; % [months]
    end

end