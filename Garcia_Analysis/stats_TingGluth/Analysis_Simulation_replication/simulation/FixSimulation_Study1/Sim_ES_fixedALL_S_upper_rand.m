% This code is a modification of Dr Chih-Chung's Code

clear all
clc
close all


% set version
version = 5

if version == 1
    m_num = 7;          % ESaDDM + z, a~OV
elseif version == 2
    m_num = 16          % aDDM, no z
elseif version == 3
    m_num = 6           % ESaDDM + z
elseif version == 4     
    m_num = 21         % ESaDDM, no z
elseif version == 5
    m_num = 2          % aDDM, z
end

% Load data
data = readtable('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv');
data.OV = data.OV_2;
data.VD = data.VD_2;
data.Choice = data.cho - 1;  %conversion
data.Correct = data.chose_right; 
cols_to_check = {'OV_2','VD_2','OV','VD','GazeDiff','FirstFixDur','FinalFixDur','MiddleFixDur', ...
                 'eachMiddleFixDur','GazeSwitch','FirstFixLoc','FinalFixLoc','DwellTimeAdvantage','corr','chose_right'};
cols_to_check = cols_to_check(ismember(cols_to_check,data.Properties.VariableNames));
data = rmmissing(data,'DataVariables',cols_to_check);

data = data(strcmp(data.phase,'ES'),:);
% if all(ismember({'DwellRight','DwellLeft'}, data.Properties.VariableNames))
%     data = data(data.DwellRight > 0 & data.DwellLeft > 0, :);
% else
%     warning('No DwellRight/DwellLeft cols');
% end
subjlist = setdiff(1:26,[1,4,5,6,14,99]);
nSubj = length(subjlist);
minSigma = randsample(0.02:0.001:0.03, nSubj, true);   % can be used

seed = sum(100*clock) + m_num + floor(1e6 * rand);
rng(seed); 

% Load posterior draws
if version == 1
    posterior_file = 'garcia_replication_For_paper_7_posterior_draws.csv';
    T = readtable(posterior_file,'VariableNamingRule','preserve');
    % simulate only first 10 draws for now
    nDraws = min(100,height(T));
    
    % Group param meta
    groupParamNames = {'a_high','a_low','a_medium','t','z', ...
                       'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S', ...
                       'theta_E','theta_S'};
    
    h = waitbar(0,'Simulating posterior draws...');
    for d = 1:nDraws
        fprintf('Draw %d/%d\n', d, nDraws);

        paramset = nan(nSubj,9);
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
    
            paraname = {['a_subj(high).' num2str(subjID)], ...
                        ['a_subj(low).' num2str(subjID)], ...
                        ['a_subj(medium).' num2str(subjID)], ...
                        ['t_subj.' num2str(subjID)], ...
                        ['z_subj.' num2str(subjID)], ...
                        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_E_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_S_subj.' num2str(subjID)]};
    
            r = zeros(1,numel(paraname));
            for k = 1:numel(paraname)
                r(k) = safeFindCol(T, paraname{k});
            end
    
            a_H = T{d,r(1)};
            a_L = T{d,r(2)};
            a_M = T{d,r(3)};
            t   = T{d,r(4)};    
            z   = T{d,r(5)};
            b2  = T{d,r(6)};
            b3  = T{d,r(7)};
            b4  = T{d,r(8)};
    
            paramset(k_subj,:) = [a_H a_L a_M t z b2 b3 b4 minSigma(k_subj)];
        end
    
        paraname_group = {'a(high)','a(low)','a(medium)','t','z', ...
                          'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S'};
        r_group = zeros(1,numel(paraname_group));
        for k = 1:numel(paraname_group)
            r_group(k) = safeFindCol(T, paraname_group{k});
        end
    
        g_aH = T{d,r_group(1)};
        g_aL = T{d,r_group(2)};
        g_aM = T{d,r_group(3)};
        g_t  = T{d,r_group(4)};
        g_z  = T{d,r_group(5)};
        g_b2 = T{d,r_group(6)};
        g_b3 = T{d,r_group(7)};
        g_b4 = T{d,r_group(8)};
    
        theta_E = g_b3 ./ g_b2;
        theta_S = g_b4 ./ g_b2;
        groupParams = [g_aH g_aL g_aM g_t g_z g_b2 g_b3 g_b4 theta_E theta_S];
    
        TBsim = [];
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
            LMH    = [1 4 5 6 7 8; 
                      2 4 5 6 7 8;
                      3 4 5 6 7 8];
            for k_OV = 1:3
                params = squeeze(paramset(k_subj,LMH(k_OV,:)));
                data_subj = data(strcmp(data.phase,'ES') & data.sub_id==subjID & data.OV==k_OV, ...
                                {'phase','sub_id','VD','p1','p2','Correct','Choice'});
                values = [data_subj.p1, data_subj.p2];
    
                behData = [];
                eyeData = [];
                for ktrial = 1:height(data_subj)
                    [Choice,RT,E,tempEyeData] = ...
                        EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);
    
                    Correct = Choice == (values(ktrial,1) < values(ktrial,2));
                    behData(ktrial,:) = [subjID ktrial k_OV data_subj.VD(ktrial) ...
                                         values(ktrial,1), values(ktrial,2) ...
                                         values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
                    eyeData(ktrial,:) = [tempEyeData.Nfix,tempEyeData.FixLocFirst,tempEyeData.FixLocLast,...
                                         tempEyeData.FixLocFirstCorr,tempEyeData.FixLocLastCorr,tempEyeData.DwellDiff,...
                                         tempEyeData.FirstFixDur,tempEyeData.MiddleFixDur,tempEyeData.FinalFixDur,...
                                         tempEyeData.eachMiddleFixDur];
                end
    
                varNames = {'sub_id','trial','OV','VD','Vl','Vr','RLdiff','Choice','Correct','rt', ...
                            'Nfix','FixLocFirst','FixLocLast','FixLocFirstCorr','FixLocLastCorr','DwellDiff', ...
                            'DwellFirst','DwellMid','DwellFinal','eachDwellMiddle'};
                temp_TBsim = array2table([behData eyeData],'VariableNames',varNames);
                TBsim = [TBsim; temp_TBsim];
            end
        end
    
        % Tag draw
        TBsim.draw_id = repmat(d,height(TBsim),1);
    
        filename = sprintf(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/' ...
                            'stats_TingGluth/Analysis_Simulation_replication/simulation/FixSimulation_Study1/ppc_full_bayesian/rand_within/rand_ESaDDM_zaOV/' ...
                            'Sim_HDDM_m7_draw%03d.mat'], d);
        save(filename,'TBsim','subjlist','groupParams','groupParamNames','seed');
    
        waitbar(d/nDraws,h);
    end
    close(h)
    

elseif version == 2
    posterior_file = 'garcia_replication_For_paper_16_posterior_draws.csv';
    T = readtable(posterior_file,'VariableNamingRule','preserve');

    nDrawsTotal = min(100, height(T));  
    startDraw   = 1;                   

    
    %Group param meta
    groupParamNames = {'a','t','z', ...
                       'v_ES_AttentionW','v_ES_InattentionW', ...
                       };
    
    h = waitbar(0,'Simulating posterior draws for the orig. aDDM...');
    for d = startDraw:nDrawsTotal
        fprintf('Draw %d/%d\n', d, nDrawsTotal);
        paramset = nan(nSubj,5);
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
    
            paraname = {['a_subj.' num2str(subjID)], ...
                        ['t_subj.' num2str(subjID)], ...
                        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_subj.' num2str(subjID)]};
    
            r = zeros(1,numel(paraname));
            for k = 1:numel(paraname)
                r(k) = safeFindCol(T, paraname{k});
            end
    
            a = T{d,r(1)};
            t   = T{d,r(2)};
            b2  = T{d,r(3)};
            b3  = T{d,r(4)};
    
            paramset(k_subj,:) = [a t b2 b3 minSigma(k_subj)];
        end
    
        paraname_group = {'a','t', ...
                          'v_ES_AttentionW','v_ES_InattentionW'};
        r_group = zeros(1,numel(paraname_group));
        for k = 1:numel(paraname_group)
            r_group(k) = safeFindCol(T, paraname_group{k});
        end
    
        g_a = T{d,r_group(1)};
        g_t  = T{d,r_group(2)};
        g_b2 = T{d,r_group(3)};
        g_b3 = T{d,r_group(4)};
    
        groupParams = [g_a g_t g_b2 g_b3];
    
        TBsim = [];
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
            LMH    = [1 2 3 4; 
                      1 2 3 4;
                      1 2 3 4];
            for k_OV = 1:3
                params = squeeze(paramset(k_subj,LMH(k_OV,:)));
                data_subj = data(strcmp(data.phase,'ES') & data.sub_id==subjID & data.OV==k_OV, ...
                                {'phase','sub_id','VD','p1','p2','Correct','Choice'});
                values = [data_subj.p1, data_subj.p2];
    
                behData = [];
                eyeData = [];
                for ktrial = 1:height(data_subj)
                    [Choice,RT,E,tempEyeData] = ...
                        EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);
    
                    Correct = Choice == (values(ktrial,1) < values(ktrial,2));
                    behData(ktrial,:) = [subjID ktrial k_OV data_subj.VD(ktrial) ...
                                         values(ktrial,1), values(ktrial,2) ...
                                         values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
                    eyeData(ktrial,:) = [tempEyeData.Nfix,tempEyeData.FixLocFirst,tempEyeData.FixLocLast,...
                                         tempEyeData.FixLocFirstCorr,tempEyeData.FixLocLastCorr,tempEyeData.DwellDiff,...
                                         tempEyeData.FirstFixDur,tempEyeData.MiddleFixDur,tempEyeData.FinalFixDur,...
                                         tempEyeData.eachMiddleFixDur];
                end
    
                varNames = {'sub_id','trial','OV','VD','Vl','Vr','RLdiff','Choice','Correct','rt', ...
                            'Nfix','FixLocFirst','FixLocLast','FixLocFirstCorr','FixLocLastCorr','DwellDiff', ...
                            'DwellFirst','DwellMid','DwellFinal','eachDwellMiddle'};
                temp_TBsim = array2table([behData eyeData],'VariableNames',varNames);
                TBsim = [TBsim; temp_TBsim];
            end
        end
    
%       tag draw
        TBsim.draw_id = repmat(d,height(TBsim),1);
    
        filename = sprintf(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/' ...
                            'stats_TingGluth/Analysis_Simulation_replication/simulation/FixSimulation_Study1/ppc_full_bayesian/rand_aDDM_noz/' ...
                            'Sim_HDDM_m16_draw%03d.mat'], d);
        save(filename,'TBsim','subjlist','groupParams','groupParamNames');
        
        progress = (d - startDraw + 1) / (nDrawsTotal - startDraw + 1);

        if ~ishandle(h)  
            h = waitbar(progress,'Simulating posterior draws for the orig. aDDM...');
        else
            waitbar(progress, h);
        end
    end

    if ishandle(h)
        close(h)
    end


elseif version == 3
    posterior_file = 'garcia_replication_For_paper_6_posterior_draws.csv';
    T = readtable(posterior_file,'VariableNamingRule','preserve');
    % simulate only first 10 draws for now
    nDraws = min(100,height(T));
    
    % Group param meta
    groupParamNames = {'a','t','z', ...
                       'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S', ...
                       'theta_E','theta_S'};
    

    h = waitbar(0,'Simulating posterior draws...');

    for d = 1:nDraws
        fprintf('Draw %d/%d\n', d, nDraws);
        paramset = nan(nSubj,7);
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
    
            paraname = {['a_subj.' num2str(subjID)], ...
                        ['t_subj.' num2str(subjID)], ...
                        ['z_subj.' num2str(subjID)], ...
                        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_E_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_S_subj.' num2str(subjID)]};
    
            r = zeros(1,numel(paraname));
            for k = 1:numel(paraname)
                r(k) = safeFindCol(T, paraname{k});
            end
    
            a   = T{d,r(1)};
            t   = T{d,r(2)};
            z   = T{d,r(3)};
            b2  = T{d,r(4)};
            b3  = T{d,r(5)};
            b4  = T{d,r(6)};
    
            paramset(k_subj,:) = [a t z b2 b3 b4 minSigma(k_subj)];
        end
    
        paraname_group = {'a','t','z', ...
                          'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S'};
        r_group = zeros(1,numel(paraname_group));
        for k = 1:numel(paraname_group)
            r_group(k) = safeFindCol(T, paraname_group{k});
        end
    
        g_a = T{d,r_group(1)};
        g_t  = T{d,r_group(2)};
        g_z  = T{d,r_group(3)};
        g_b2 = T{d,r_group(4)};
        g_b3 = T{d,r_group(5)};
        g_b4 = T{d,r_group(6)};
    
        theta_E = g_b3 ./ g_b2;
        theta_S = g_b4 ./ g_b2;
        groupParams = [g_a g_t g_z g_b2 g_b3 g_b4 theta_E theta_S];
    
        TBsim = [];
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
            LMH    = [1 2 3 4 5 6; 
                      1 2 3 4 5 6;
                      1 2 3 4 5 6];
            for k_OV = 1:3
                params = squeeze(paramset(k_subj,LMH(k_OV,:)));
                data_subj = data(strcmp(data.phase,'ES') & data.sub_id==subjID & data.OV==k_OV, ...
                                {'phase','sub_id','VD','p1','p2','Correct','Choice'});
                values = [data_subj.p1, data_subj.p2];
    
                behData = [];
                eyeData = [];
                for ktrial = 1:height(data_subj)
                    [Choice,RT,E,tempEyeData] = ...
                        EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);
    
                    Correct = Choice == (values(ktrial,1) < values(ktrial,2));
                    behData(ktrial,:) = [subjID ktrial k_OV data_subj.VD(ktrial) ...
                                         values(ktrial,1), values(ktrial,2) ...
                                         values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
                    eyeData(ktrial,:) = [tempEyeData.Nfix,tempEyeData.FixLocFirst,tempEyeData.FixLocLast,...
                                         tempEyeData.FixLocFirstCorr,tempEyeData.FixLocLastCorr,tempEyeData.DwellDiff,...
                                         tempEyeData.FirstFixDur,tempEyeData.MiddleFixDur,tempEyeData.FinalFixDur,...
                                         tempEyeData.eachMiddleFixDur];
                end
    
                varNames = {'sub_id','trial','OV','VD','Vl','Vr','RLdiff','Choice','Correct','rt', ...
                            'Nfix','FixLocFirst','FixLocLast','FixLocFirstCorr','FixLocLastCorr','DwellDiff', ...
                            'DwellFirst','DwellMid','DwellFinal','eachDwellMiddle'};
                temp_TBsim = array2table([behData eyeData],'VariableNames',varNames);
                TBsim = [TBsim; temp_TBsim];
            end
        end
    
        % Tag draw
        TBsim.draw_id = repmat(d,height(TBsim),1);
    
        filename = sprintf(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/' ...
                            'stats_TingGluth/Analysis_Simulation_replication/simulation/FixSimulation_Study1/ppc_full_bayesian/rand_within/rand_ESaDDM_z/' ...
                            'Sim_HDDM_m6_draw%03d.mat'], d);
        save(filename,'TBsim','subjlist','groupParams','groupParamNames');
    
        waitbar(d/nDraws,h);
    end
    close(h)

elseif version == 4
    posterior_file = 'garcia_replication_For_paper_21_posterior_draws.csv';
    T = readtable(posterior_file,'VariableNamingRule','preserve');
    % simulate only first 10 draws for now
    nDraws = min(100,height(T));
    
    groupParamNames = {'a','t', ...
                       'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S', ...
                       'theta_E','theta_S'};
    
    % loop and save draw independently
    h = waitbar(0,'Simulating posterior draws...');
    for d = 1:nDraws
        fprintf('Draw %d/%d\n', d, nDraws);
        paramset = nan(nSubj,6);
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
    
            paraname = {['a_subj.' num2str(subjID)], ...
                        ['t_subj.' num2str(subjID)], ...
                        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_E_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_S_subj.' num2str(subjID)]};
    
            r = zeros(1,numel(paraname));
            for k = 1:numel(paraname)
                r(k) = safeFindCol(T, paraname{k});
            end
    
            a   = T{d,r(1)};
            t   = T{d,r(2)};
            b2  = T{d,r(3)};
            b3  = T{d,r(4)};
            b4  = T{d,r(5)};
    
            paramset(k_subj,:) = [a t b2 b3 b4 minSigma(k_subj)];
        end
    
        paraname_group = {'a','t', ...
                          'v_ES_AttentionW','v_ES_InattentionW_E','v_ES_InattentionW_S'};
        r_group = zeros(1,numel(paraname_group));
        for k = 1:numel(paraname_group)
            r_group(k) = safeFindCol(T, paraname_group{k});
        end
    
        g_a = T{d,r_group(1)};
        g_t  = T{d,r_group(2)};
        g_b2 = T{d,r_group(3)};
        g_b3 = T{d,r_group(4)};
        g_b4 = T{d,r_group(5)};
    
        theta_E = g_b3 ./ g_b2;
        theta_S = g_b4 ./ g_b2;
        groupParams = [g_a g_t g_b2 g_b3 g_b4 theta_E theta_S];
    
        TBsim = [];
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
            LMH    = [1 2 3 4 5; 
                      1 2 3 4 5;
                      1 2 3 4 5];
            for k_OV = 1:3
                params = squeeze(paramset(k_subj,LMH(k_OV,:)));
                data_subj = data(strcmp(data.phase,'ES') & data.sub_id==subjID & data.OV==k_OV, ...
                                {'phase','sub_id','VD','p1','p2','Correct','Choice'});
                values = [data_subj.p1, data_subj.p2];
    
                behData = [];
                eyeData = [];
                for ktrial = 1:height(data_subj)
                    [Choice,RT,E,tempEyeData] = ...
                        EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);
    
                    Correct = Choice == (values(ktrial,1) < values(ktrial,2));
                    behData(ktrial,:) = [subjID ktrial k_OV data_subj.VD(ktrial) ...
                                         values(ktrial,1), values(ktrial,2) ...
                                         values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
                    eyeData(ktrial,:) = [tempEyeData.Nfix,tempEyeData.FixLocFirst,tempEyeData.FixLocLast,...
                                         tempEyeData.FixLocFirstCorr,tempEyeData.FixLocLastCorr,tempEyeData.DwellDiff,...
                                         tempEyeData.FirstFixDur,tempEyeData.MiddleFixDur,tempEyeData.FinalFixDur,...
                                         tempEyeData.eachMiddleFixDur];
                end
    
                varNames = {'sub_id','trial','OV','VD','Vl','Vr','RLdiff','Choice','Correct','rt', ...
                            'Nfix','FixLocFirst','FixLocLast','FixLocFirstCorr','FixLocLastCorr','DwellDiff', ...
                            'DwellFirst','DwellMid','DwellFinal','eachDwellMiddle'};
                temp_TBsim = array2table([behData eyeData],'VariableNames',varNames);
                TBsim = [TBsim; temp_TBsim];
            end
        end
    
        % Tag draw
        TBsim.draw_id = repmat(d,height(TBsim),1);
    
        filename = sprintf(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/' ...
                            'stats_TingGluth/Analysis_Simulation_replication/simulation/FixSimulation_Study1/ppc_full_bayesian/rand_ESaDDM_noz/' ...
                            'Sim_HDDM_m21_draw%03d.mat'], d);
        save(filename,'TBsim','subjlist','groupParams','groupParamNames');
    
        waitbar(d/nDraws,h);
    end
    close(h)

elseif version == 5
   posterior_file = 'garcia_replication_For_paper_2_posterior_draws.csv';
    T = readtable(posterior_file,'VariableNamingRule','preserve');
    % simulate only first 10 draws for now
    nDraws = min(100,height(T));
    
    % Group param meta
    groupParamNames = {'a','t','z', ...
                       'v_ES_AttentionW','v_ES_InattentionW', ...
                       'theta'};
    
    h = waitbar(0,'Simulating posterior draws...');

    for d = 1:nDraws
        fprintf('Draw %d/%d\n', d, nDraws);
        paramset = nan(nSubj,6);
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
    
            paraname = {['a_subj.' num2str(subjID)], ...
                        ['t_subj.' num2str(subjID)], ...
                        ['z_subj.' num2str(subjID)], ...
                        ['v_ES_AttentionW_subj.' num2str(subjID)], ...
                        ['v_ES_InattentionW_subj.' num2str(subjID)]};
    
            r = zeros(1,numel(paraname));
            for k = 1:numel(paraname)
                r(k) = safeFindCol(T, paraname{k});
            end
    
            a   = T{d,r(1)};
            t   = T{d,r(2)};
            z   = T{d,r(3)};
            b2  = T{d,r(4)};
            b3  = T{d,r(5)};
    
            paramset(k_subj,:) = [a t z b2 b3 minSigma(k_subj)];
        end
    
        paraname_group = {'a','t','z', ...
                          'v_ES_AttentionW','v_ES_InattentionW'};
        r_group = zeros(1,numel(paraname_group));
        for k = 1:numel(paraname_group)
            r_group(k) = safeFindCol(T, paraname_group{k});
        end
    
        g_a = T{d,r_group(1)};
        g_t  = T{d,r_group(2)};
        g_z  = T{d,r_group(3)};
        g_b2 = T{d,r_group(4)};
        g_b3 = T{d,r_group(5)};
    
        theta= g_b3 ./ g_b2;
        groupParams = [g_a g_t g_z g_b2 g_b3 theta];
    
        TBsim = [];
        for k_subj = 1:nSubj
            subjID = subjlist(k_subj);
            LMH    = [1 2 3 4 5; 
                      1 2 3 4 5;
                      1 2 3 4 5];
            for k_OV = 1:3
                params = squeeze(paramset(k_subj,LMH(k_OV,:)));
                data_subj = data(strcmp(data.phase,'ES') & data.sub_id==subjID & data.OV==k_OV, ...
                                {'phase','sub_id','VD','p1','p2','Correct','Choice'});
                values = [data_subj.p1, data_subj.p2];
    
                behData = [];
                eyeData = [];
                for ktrial = 1:height(data_subj)
                    [Choice,RT,E,tempEyeData] = ...
                        EvidenceAccumulate_S_upper(k_OV, values(ktrial,1), values(ktrial,2), params);
    
                    Correct = Choice == (values(ktrial,1) < values(ktrial,2));
                    behData(ktrial,:) = [subjID ktrial k_OV data_subj.VD(ktrial) ...
                                         values(ktrial,1), values(ktrial,2) ...
                                         values(ktrial,2)-values(ktrial,1) Choice Correct RT/1000];
                    eyeData(ktrial,:) = [tempEyeData.Nfix,tempEyeData.FixLocFirst,tempEyeData.FixLocLast,...
                                         tempEyeData.FixLocFirstCorr,tempEyeData.FixLocLastCorr,tempEyeData.DwellDiff,...
                                         tempEyeData.FirstFixDur,tempEyeData.MiddleFixDur,tempEyeData.FinalFixDur,...
                                         tempEyeData.eachMiddleFixDur];
                end
    
                varNames = {'sub_id','trial','OV','VD','Vl','Vr','RLdiff','Choice','Correct','rt', ...
                            'Nfix','FixLocFirst','FixLocLast','FixLocFirstCorr','FixLocLastCorr','DwellDiff', ...
                            'DwellFirst','DwellMid','DwellFinal','eachDwellMiddle'};
                temp_TBsim = array2table([behData eyeData],'VariableNames',varNames);
                TBsim = [TBsim; temp_TBsim];
            end
        end
    
        % Tag draw
        TBsim.draw_id = repmat(d,height(TBsim),1);
    
        filename = sprintf(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/' ...
                            'stats_TingGluth/Analysis_Simulation_replication/simulation/FixSimulation_Study1/ppc_full_bayesian/rand_aDDM_z/' ...
                            'Sim_HDDM_m2_draw%03d.mat'], d);
        save(filename,'TBsim','subjlist','groupParams','groupParamNames');
    
        waitbar(d/nDraws,h);
    end
    close(h)
end

% helper
% ---------------------------------------------------------------------------------

function colIdx = safeFindCol(T, expectedName)
    % variable names from table
    varNames = T.Properties.VariableNames;
    expectedNorm = regexprep(expectedName, '[^a-zA-Z0-9_]', '');
    colIdx = [];
    for i = 1:length(varNames)
        thisName = regexprep(varNames{i}, '[^a-zA-Z0-9_]', '');
        % check if it matches
        if strcmp(thisName, expectedNorm)
            colIdx = i;
            break
        end
    end
    if isempty(colIdx)
        error(['Column ', expectedName, ' not found']);
    end
end

