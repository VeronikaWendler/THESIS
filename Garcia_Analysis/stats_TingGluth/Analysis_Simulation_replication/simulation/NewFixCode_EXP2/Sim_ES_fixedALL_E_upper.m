% IMPORTANT: This code assumes response coding with the upper bound being
% the E option (for analysisi purposes E is always the left option; so when
% chose_left == 1, E was chosen
% Modified code from Chih-Chung Ting

clear all
clc
close all

m_num = 28
seed = sum(100*clock) + m_num + floor(1e6 * rand);
rng(seed);

%% Load task-related information
%-----------------------------------------

data = readtable('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv');
                 
data.OV = data.OV_2;
data.VD = data.VD_2;

% Convert 'cho' 1 = left, 2 = right into 'Choice' 1 = left, 0 = right
data.Choice = abs(data.cho - 2);
data.Choice

% Rename 'corr' to 'Correct'
% for accuracy coding corr ==1 if correct, 0 otherwise
% chose_left is 1 if left was chosen (E option chosen), 0 otherwise (S
% option) used for response coding - chose_left is like our response col
% here
data.Correct = data.chose_left;

cols_to_check = {'OV_2', 'VD_2', 'OV', 'VD', 'GazeDiff', 'FirstFixDur', 'FinalFixDur', 'MiddleFixDur', ...
                 'eachMiddleFixDur', 'GazeSwitch', 'FirstFixLoc', 'FinalFixLoc', ...
                 'DwellTimeAdvantage', 'chose_right'};

existing_vars = ismember(cols_to_check, data.Properties.VariableNames);
cols_to_check = cols_to_check(existing_vars); 

data.OV = data.OV_2;
data.VD = data.VD_2;

data = rmmissing(data, 'DataVariables', cols_to_check);

data = data(strcmp(data.phase, 'ES'), :);

subjlist = setdiff([1:26], [6, 14, 20, 26, 2, 9, 18]);
nSubj = length(subjlist);
minSigma = randsample(0.02:0.001:0.03, nSubj, true); 

data = data(strcmp(data.phase, 'ES'), :);

m_num = 28; 
param_file = ['params_ES_VAL_m' num2str(m_num) '.csv'];
converging_file = ['gelman_rubin_ES_VAL_m' num2str(m_num) '.csv'];
conv_criteria = 1.1;

%% Load parameter sets
T = readtable(fullfile(param_file));
Tconv = readtable(fullfile(converging_file));

%% load individual parameter set
%-----------------------------------------
% for k_subj = 1:nSubj
%     clear paraname r
%     subjID = subjlist(k_subj);
%     paraname = {['a_subj.' num2str(subjID)], ...
%         ['t_subj(high).' num2str(subjID)], ...
%         ['t_subj(low).' num2str(subjID)], ...
%         ['t_subj(medium).' num2str(subjID)], ...
%         ['z_subj.' num2str(subjID)], ...
%         ['v_Intercept_subj.' num2str(subjID)], ...
%         ['v_ES_AttentionW_subj.' num2str(subjID)], ...
%         ['v_ES_InattentionW:C(OVcate)[high]_subj.' num2str(subjID)], ...
%         ['v_ES_InattentionW:C(OVcate)[low]_subj.' num2str(subjID)], ...
%         ['v_ES_InattentionW:C(OVcate)[medium]_subj.' num2str(subjID)], ...
%         };
%     
%     r = find(ismember(T.Var1, paraname));
%     a =   T.mean(r(1));
%     t_H = T.mean(r(2));
%     t_L = T.mean(r(3));
%     t_M = T.mean(r(4));
%     z = T.mean(r(5));
%     beta0 = T.mean(r(6));
%     beta2 = T.mean(r(7));
%     beta3_H = T.mean(r(8));
%     beta3_L = T.mean(r(9));
%     beta3_M = T.mean(r(10));
% 
%     paramset(k_subj,:) = [a t_H t_L t_M z beta0 beta2 beta3_H beta3_L beta3_M minSigma(k_subj)];
%     
%         paraname_group_name = {'a', ...
%                 't(high)', ...
%                 't(low)', ...
%                 't(medium)', ...
%                 'z', ...
%                 'v_Intercept', ...
%                 'v_ES_AttentionW', ...
%                 'v_ES_InattentionW:C(OVcate)[high]', ...
%                 'v_ES_InattentionW:C(OVcate)[low]', ...
%                 'v_ES_InattentionW:C(OVcate)[medium]'
%         };
% 
%     r_group = find(ismember(T.Var1, paraname_group_name));
%     a =   T.mean(r_group(1));
%     t_H = T.mean(r_group(2));
%     t_L = T.mean(r_group(3));
%     t_M = T.mean(r_group(4));
%     z = T.mean(r_group(5));
%     beta0 = T.mean(r_group(6));
%     beta2 = T.mean(r_group(7));
%     beta3_H = T.mean(r_group(8));
%     beta3_L = T.mean(r_group(9));
%     beta3_M = T.mean(r_group(10));
% 
%     theta(1,1) = beta3_H./beta2;
%     theta(1,2) = beta3_L./beta2;
%     theta(1,3) = beta3_M./beta2;
% 
%     paraname_group = [a t_H t_L t_M z beta0 beta2 beta3_H beta3_L beta3_M theta];
% end
% convergeset = Tconv.Gelman_Rubin;


% for k_subj = 1:nSubj
%     subjID = subjlist(k_subj);
%     paraname = {['a_subj.' num2str(subjID)], ...
%                 ['t_subj.' num2str(subjID)], ...
%                 ['z_subj.' num2str(subjID)], ...
%                 ['v_Intercept_subj.' num2str(subjID)], ...
%                 ['v_AttentionW_E_subj.' num2str(subjID)], ...
%                 ['v_AttentionW_S_subj.' num2str(subjID)], ...
%                 ['v_InattentionW_E_subj.' num2str(subjID)], ...
%                 ['v_InattentionW_S_subj.' num2str(subjID)]};
%     r = find(ismember(T.Var1, paraname));
% 
%     a     = T.mean(r(1));
%     ndt   = T.mean(r(2));
%     z     = T.mean(r(3));
%     beta0 = T.mean(r(4));
%     beta1 = T.mean(r(5));
%     beta2 = T.mean(r(6));
%     beta3 = T.mean(r(7));
%     beta4 = T.mean(r(8));
% 
%     paramset(k_subj,:) = [a ndt z beta0 beta1 beta2 beta3 beta4 minSigma(k_subj)];
% 
%     % group-level (optional; you saved it before, keep parity)
%     paraname_group_name = {'a', ...
%         't', ...
%         'z', ...
%         'v_Intercept', ...
%         'v_AttentionW_E', ...
%         'v_AttentionW_S', ...
%         'v_InattentionW_E', ...
%         'v_InattentionW_S'}; ...
% 
%     r_group = find(ismember(T.Var1, paraname_group_name));
%     a     = T.mean(r_group(1));
%     ndt   = T.mean(r_group(2));
%     z   = T.mean(r_group(3));
%     beta0 = T.mean(r_group(4));
%     beta1 = T.mean(r_group(5));
%     beta2 = T.mean(r_group(6));
%     beta3 = T.mean(r_group(7));
%     beta4 = T.mean(r_group(8));
%     theta_E(1,1) = beta3./beta1;
%     theta_S(1,1) = beta4./beta2;
%     paraname_group = [a ndt z beta0 beta1 beta2 beta3 beta4 theta_E theta_S];
% end
%     convergeset = Tconv.Gelman_Rubin;

for k_subj = 1:nSubj
    clear paraname r
    subjID = subjlist(k_subj);
    paraname = {['a_subj.' num2str(subjID)], ...
        ['t_subj.' num2str(subjID)], ...
        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
        ['v_ES_InattentionW_subj.' num2str(subjID)], ...
        };
    
    r = find(ismember(T.Var1, paraname));
    a =   T.mean(r(1));
    t = T.mean(r(2));
    beta0 = T.mean(r(3));
    beta1 = T.mean(r(4));

    paramset(k_subj,:) = [a t beta0 beta1 minSigma(k_subj)];
    
        paraname_group_name = {'a', ...
                't', ...
                'v_ES_AttentionW', ...
                'v_ES_InattentionW', ...
        };

    r_group = find(ismember(T.Var1, paraname_group_name));
    a =   T.mean(r_group(1));
    t = T.mean(r_group(2));
    beta0 = T.mean(r_group(3));
    beta1 = T.mean(r_group(4));

    theta(1,1) = beta1./beta0;
   
    paraname_group = [a t beta0 beta1 theta];
end
convergeset = Tconv.Gelman_Rubin;

save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode/Sim/Sim_HDDM_ES_VAL_m28_newfix_corr_accumfun_2_ceqc','paramset','convergeset','paraname_group')

% Nrep = 100;

% Start simulation
load('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode/Sim/Sim_HDDM_ES_VAL_m28_newfix_corr_accumfun_2_ceqc')
TBsim = [];
h = waitbar(0, 'Please wait...');
for k_subj = 1:nSubj
    count = 0;
    temp_TBsim = [];
    subjID = subjlist(k_subj);
    LMH    = [1 2 3 4 ; ...
        1 2 3 4 ; ...
        1 2 3 4 ];

    waitbar(k_subj / nSubj)
    
    for k_OV = 1:3
        display(['subj' num2str(subjID) '_OV' num2str(k_OV)])
        params = squeeze(paramset(k_subj,LMH(k_OV,:)));
        
        data.OV = data.OV_2;
        data.VD = data.VD_2;

        % Filter data for this subject & OV level
        data_subj = data(strcmp(data.phase, 'ES') & data.sub_id == subjID & data.OV == k_OV, {'phase', 'sub_id', 'VD', 'p1', 'p2', 'Correct', 'Choice'});
        values = [data_subj.p1, data_subj.p2].*1;


        for ktrial = 1:length(data_subj.sub_id)
            % Simulate aDDM process
            [Choice, RT, E, tempEyeData, FixAaLL, sumdvALL, FixDur] = EvidenceAccumulate_E_upper(k_OV, values(ktrial,1), values(ktrial,2), params);

            Correct = Choice == (values(ktrial,1) > values(ktrial,2));   
            
            count = count + 1;
            behData(count,:) = [subjID+1000 ktrial k_OV data_subj.VD(ktrial) values(ktrial,1),values(ktrial,2) values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
            eyeData(count,:) = [tempEyeData.Nfix, tempEyeData.FixLocFirst, tempEyeData.FixLocLast,tempEyeData.FixLocFirstCorr, tempEyeData.FixLocLastCorr, tempEyeData.DwellDiff, ...
                tempEyeData.FirstFixDur, tempEyeData.MiddleFixDur, tempEyeData.FinalFixDur, tempEyeData.eachMiddleFixDur];
        end
        
        varNames = {'sub_id';'trial';'OV';'VD';'Vl';'Vr';'RLdiff';'Choice';'Correct';'rt'; ...
            'Nfix';'FixLocFirst';'FixLocLast'; 'FixLocFirstCorr';'FixLocLastCorr'; 'DwellDiff'; ...
            'DwellFirst';'DwellMid';'DwellFinal';'eachDwellMiddle'};
        temp_TBsim = table(behData(:,1),behData(:,2),behData(:,3),behData(:,4),behData(:,5),behData(:,6), behData(:,7),behData(:,8),behData(:,9),behData(:,10), ...
            eyeData(:,1),eyeData(:,2),eyeData(:,3),eyeData(:,4),eyeData(:,5),eyeData(:,6), ...
            eyeData(:,7),eyeData(:,8),eyeData(:,9),eyeData(:,10),'VariableNames',varNames);
    end
    TBsim = [TBsim;temp_TBsim];
end

save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode/Sim/sim_HDDM_ES_VAL_m28_newfix_corr_accumfun_2_ceqc','TBsim','paramset','subjlist')
close(h)




