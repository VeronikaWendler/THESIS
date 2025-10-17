%% Clean up the environment
%-----------------------------------------
% IMPORTANT: This code assumes response coding with the upper bound being
% the S option (for analysisi purposes E is always the left option; so when
% chose_left == 1, E was chosen, when chose_right == 1, S was chosen

clear all
clc
close all

m_num = 6
seed = sum(100*clock) + m_num + floor(1e6 * rand);
rng(seed);


data = readtable('D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv');
                 
data.OV = data.OV_2;
data.VD = data.VD_2;

% Convert 'cho' 1 = left, 2 = right into 'Choice' 1 = right, 0 = left
data.Choice = data.cho - 1;
data.Choice

% Rename 'corr' to 'Correct'
% for accuracy coding corr ==1 if correct, 0 otherwise
% chose_left is 1 if left was chosen (E option chosen), 0 otherwise (S
% option) used for response coding - chose_left is like our response col
% here
data.Correct = data.chose_right;

cols_to_check = {'OV_2', 'VD_2', 'OV', 'VD', 'GazeDiff', 'FirstFixDur', 'FinalFixDur', 'MiddleFixDur', ...
                 'eachMiddleFixDur', 'GazeSwitch', 'FirstFixLoc', 'FinalFixLoc', ...
                 'DwellTimeAdvantage', 'chose_right'};

existing_vars = ismember(cols_to_check, data.Properties.VariableNames);
cols_to_check = cols_to_check(existing_vars); 

data.OV = data.OV_2;
data.VD = data.VD_2;

data = rmmissing(data, 'DataVariables', cols_to_check);
data = data(strcmp(data.phase, 'ES'), :);
subjlist = setdiff([1:27], [6, 14, 20, 26, 2, 9, 18]);
nSubj = length(subjlist);
minSigma = randsample(0.02:0.001:0.03, nSubj, true); % HDDM assumes intra-trial variance is 1 for each time unit.
data = data(strcmp(data.phase, 'ES'), :);

m_num = 6; 
param_file = ['params_For_paper_m' num2str(m_num) '.csv'];
converging_file = ['gelman_rubin_For_paper_m' num2str(m_num) '.csv'];
conv_criteria = 1.1;


%% Load parameter sets
T = readtable(fullfile(param_file));
Tconv = readtable(fullfile(converging_file));

%% load individual parameter set

for k_subj = 1:nSubj
    clear paraname r
    subjID = subjlist(k_subj);
    paraname = {['a_subj(high).' num2str(subjID)], ...
        ['a_subj(low).' num2str(subjID)], ...
        ['a_subj(medium).' num2str(subjID)], ...
        ['t_subj.' num2str(subjID)], ...
        ['z_subj.' num2str(subjID)], ...
        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
        ['v_ES_InattentionW_E_subj.' num2str(subjID)], ...
        ['v_ES_InattentionW_S_subj.' num2str(subjID)], ...
        };
    
    r = find(ismember(T.Var1, paraname));
    a_H =   T.mean(r(1));
    a_L =   T.mean(r(2));
    a_M =   T.mean(r(3));
    t =     T.mean(r(4));
    z =     T.mean(r(5));
    beta2 = T.mean(r(6));
    beta3 = T.mean(r(7));
    beta4 = T.mean(r(8));


    paramset(k_subj,:) = [a_H a_L a_M t z beta2 beta3 beta4 minSigma(k_subj)];
    
        paraname_group_name = {'a(high)', ...
                'a(low)', ...
                'a(medium)', ...
                't', ...
                'z', ...
                'v_ES_AttentionW', ...
                'v_ES_InattentionW_E', ...
                'v_ES_InattentionW_S', ...
        };

    r_group = find(ismember(T.Var1, paraname_group_name));
    a_H =   T.mean(r_group(1));
    a_L =   T.mean(r_group(2));
    a_M =   T.mean(r_group(3));
    t = T.mean(r_group(4));
    z = T.mean(r_group(5));
    beta2 = T.mean(r_group(6));
    beta3 = T.mean(r_group(7));
    beta4 = T.mean(r_group(8));


%     theta = beta3./beta2;
     theta_E = beta3./beta2;
     theta_S = beta4./beta2;

    paraname_group = [a_H a_L a_M t z beta2 beta3 beta4 theta_E theta_S];


end
convergeset = Tconv.Gelman_Rubin;

save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/Sim_ESaDDM/Sim_HDDM_For_paper_m6_sim10','paramset','convergeset','paraname_group')

% Nrep = 100;

%% Start simulation
load('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/Sim_ESaDDM/Sim_HDDM_For_paper_m6_sim10')
TBsim = [];
h = waitbar(0, 'Please wait...');
for k_subj = 1:nSubj
    count = 0;
    temp_TBsim = [];
    subjID = subjlist(k_subj);
    LMH    = [1 4 5 6 7 8; ...
        2 4 5 6 7 8; ...
        3 4 5 6 7 8];

    waitbar(k_subj / nSubj)
    
    for k_OV = 1:3
        display(['subj' num2str(subjID) '_OV' num2str(k_OV)])
        params = squeeze(paramset(k_subj,LMH(k_OV,:)));
        
        data.OV = data.OV_2;
        data.VD = data.VD_2;

        % Filter data for this subject & OV level
        data_subj = data(strcmp(data.phase, 'ES') & data.sub_id == subjID & data.OV == k_OV, {'phase', 'sub_id', 'VD', 'p1', 'p2', 'Correct', 'Choice'});
        values = [data_subj.p1, data_subj.p2].*1; % Check this ... Convert values to percentages (not needed, only if your values are decimals do *100)


        for ktrial = 1:length(data_subj.sub_id)
            % Simulate aDDM process
            [Choice, RT, E, tempEyeData, FixAaLL, sumdvALL, FixDur] = EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);

            Correct = Choice == (values(ktrial,1) < values(ktrial,2));    % == (values(ktrial,1) < values(ktrial,2))                     == (values(ktrial,1) > values(ktrial,2)) including this is labeled as 2_ceqc
            
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

save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/Sim_ESaDDM/sim_HDDM_For_paper_m6_sim10','TBsim','paramset','subjlist')
close(h)
