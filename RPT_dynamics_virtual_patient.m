function dYdt = RPT_dynamics_virtual_patient(t, Y, param_current_VP, Ee_Lu177)
%--------------------------------------------------------------------------
% RPT_dynamics_virtual_patient
%
% Author: Matteo Italia, ORCID ID 0000-0001-9633-8396 (email: matteo.italia@uclm.es)
%
% Defines the system of ODEs describing tumor growth and radiopharmaceutical
% kinetics in multiple organs during 177Lu-PSMA therapy between consecutive activity injections.
%
% The model accounts for tumor proliferation, sublethal and lethal damage,
% and activity clearance from tumor and normal tissues.
%
%
% INPUTS:
%   t                 - time (h)
%   Y                 - state vector:
%                       [MT_P, MT_SD, MT_D, A_tumor, A_kidney,
%                        A_parotid, A_lacrimal, A_whole_body, A_marrow, AD_tumor]
%   param_current_VP - vector of virtual patient parameters
%
% OUTPUT:
%   dYdt              - derivatives of Y at time t
%
%--------------------------------------------------------------------------
%% Extract patient-specific parameters
M_TP0          = param_current_VP(1, 1);
MK             = param_current_VP(1, 2);
MPG            = param_current_VP(1, 3);
MLG            = param_current_VP(1, 4);
MRM            = param_current_VP(1, 5);
alpha          = param_current_VP(1, 6);
beta           = param_current_VP(1, 7);
p              = param_current_VP(1, 8);
Epsilon        = param_current_VP(1, 9);
muT            = param_current_VP(1, 10);
kc             = param_current_VP(1, 11);
kg             = param_current_VP(1, 12);
lambda_phy     = param_current_VP(1, 13);
lambda_r_T     = param_current_VP(1, 14);
lambda_r_K     = param_current_VP(1, 15);
lambda_r_PG    = param_current_VP(1, 16);
lambda_r_LG    = param_current_VP(1, 17);
lambda_r_RM    = param_current_VP(1, 18);
lambda_r_RoB    = param_current_VP(1, 19);
alpha_over_beta_K = param_current_VP(1, 20);
mu_K           = param_current_VP(1, 21);

%% Unpack state variables
M_TP        = Y(1);  % Proliferating tumor mass (g)
M_TSD       = Y(2);  % Sublethally damaged tumor mass (g)
M_TD        = Y(3);  % Dead tumor mass (g)
A_tumor     = Y(4);  % Tumor activity (GBq)
A_kidney    = Y(5);  % Kidney activity (GBq)
A_parotid   = Y(6);  % Parotid gland activity (GBq)
A_lacrimal  = Y(7);  % Lacrimal gland activity (GBq)
A_rest_body= Y(8);  % Rest of the body activity (GBq)
A_marrow    = Y(9);  % Red marrow activity (GBq)
AD_tumor    = Y(10); % Tumor absorbed dose (Gy)

%% ------------------------------------------------------------------------
%                           TUMOR DYNAMICS
% -------------------------------------------------------------------------
% Dose rate to tumor 
D_tumor = Ee_Lu177 * (M_TP + M_TSD + M_TD)^(-0.99) * A_tumor;

% Damage rate constants (h⁻¹)
k_alpha = alpha * D_tumor;               % Lethal damage rate
k_p     = 2 * p * D_tumor;               % Sublethal damage rate
k_pp    = (p + alpha) * Epsilon * D_tumor; % Sublethal interaction leading to lethally damaged

% ODEs for tumor compartments
dM_TPdt  = kg * M_TP - (k_p + k_alpha) * M_TP + muT * M_TSD;
dM_TSDdt = k_p * M_TP - (muT + k_pp) * M_TSD;
dM_TDdt  = k_alpha * M_TP + k_pp * M_TSD - kc * M_TD;

%% ------------------------------------------------------------------------
%                ORGAN AND TUMOR ACTIVITY KINETICS
% -------------------------------------------------------------------------
% Simple mono-exponential clearance for each compartment
dA_tumordt     = - (lambda_phy + lambda_r_T)  * A_tumor;
dA_kidneydt    = - (lambda_phy + lambda_r_K)  * A_kidney;
dA_parotiddt   = - (lambda_phy + lambda_r_PG) * A_parotid;
dA_lacrimaldt  = - (lambda_phy + lambda_r_LG) * A_lacrimal;
dA_whole_bodydt= - (lambda_phy + lambda_r_RoB) * A_rest_body;
dA_marrowdt    = - (lambda_phy + lambda_r_RM) * A_marrow;

%% ------------------------------------------------------------------------
%                     ABSORBED DOSE ACCUMULATION
% -------------------------------------------------------------------------
dAD_tumor = D_tumor; % Integrate absorbed dose over time

%% ------------------------------------------------------------------------
%                            OUTPUT VECTOR
% -------------------------------------------------------------------------
dYdt = [dM_TPdt; dM_TSDdt; dM_TDdt; ...
        dA_tumordt; dA_kidneydt; dA_parotiddt; dA_lacrimaldt; ...
        dA_whole_bodydt; dA_marrowdt; dAD_tumor];

end