% analyze simulated data with simple statistics (e.g., ttest and correlations)

clear all
clc
close all

%% list of data
DataName   = {'realData' 'simData'};   

TempRealData = load('Data/GarciaData_eye_organization_2.mat');

answer = questdlg('Data of interest?', ...
    'Plot Menu', ...
    'ES', 'ES');


switch answer
    case 'EE'
        TB.realData = TempRealData.TB.Garcia(strcmp(TempRealData.TB.Garcia.phase, 'EE'), :);
        
        TB.realData.OV = TB.realData.OV_2;
        TB.realData.VD = TB.realData.VD_2;

        excludeSubs = [1, 4, 5, 6, 14, 99];
        TB.realData = TempRealData.TB.Garcia(~ismember(TempRealData.TB.Garcia.sub_id, excludeSubs), :);

        TB.simData = [];

        TempSimData = load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Sim/Value_mod/sim_HDDM_ES_VAL_m1.mat']);
        TempSimData2= load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Sim/without_Fix_weight_FixR/sim_ES_ZBIAS_sig03.mat']);                                      %sim_Bright_fixedALL sim_Bright_fixThreashold sim_Bright_sig025 sim_Bright_fixTheta sim_Bright_variedA
        realData = 1;
        ylimSet = [0 3; 50 100]; % 
        
    case 'ES'
        TB.realData = TempRealData.TB.Garcia(strcmp(TempRealData.TB.Garcia.phase,'ES'), :);
        TB.realData.OV = TB.realData.OV_2;
        TB.realData.VD = TB.realData.VD_2;
        excludeSubs = [1,4,5,6,14,99];
        TB.realData = TB.realData(~ismember(TB.realData.sub_id, excludeSubs), :);
        
        TB.simData = [];
        TempSimData = load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode/Sim/sim_HDDM_ES_VAL_m37_newfix_corr_accumfun_2_ceqc.mat']);
        %TempSimData2= load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Sim/without_Fix_weight_FixR/sim_ES_ZBIAS_m14.mat']);                                  
        realData = 2;
        ylimSet = [0 3; 50 100]; % 

    case 'ESEE'
        TB.realData = TempRealData.TB.Garcia(ismember(TempRealData.TB.Garcia.phase, {'ES','EE'}), :);

        excludeSubs = [1, 4, 5, 6, 14, 99];
        TB.realData = TempRealData.TB.Garcia(~ismember(TempRealData.TB.Garcia.sub_id, excludeSubs), :);
        
        TB.simData = [];
        TempSimData = load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Sim/with_Fix_weight_without_FixR/sim_ESEE_fixedALL.mat']);
        TempSimData2= load(['D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Sim/with_Fix_weight_without_FixR/sim_ESEE_5_sig03.mat']);                                  
        realData = 3;
        ylimSet = [0 3; 50 100]; % 

end


TB.simData   = TempSimData.TBsim(TempSimData.TBsim.Nfix>=0,:);
%TB.simData2  = TempSimData2.TBsim(TempSimData2.TBsim.Nfix>=0,:);
nSubj = [length(unique(TB.realData.sub_id)) length(unique(TB.simData.sub_id)) ];   %length(unique(TB.simData2.sub_id))


if strcmp(answer,'ESEE')
    if ismember('phase', TB.simData.Properties.VariableNames) && isnumeric(TB.simData.phase)
        TB.simData.phase = string(TB.simData.phase);
        TB.simData.phase(TB.simData.phase=="1") = "ES";
        TB.simData.phase(TB.simData.phase=="2") = "EE";
    end
%     if ismember('phase', TB.simData2.Properties.VariableNames) && isnumeric(TB.simData2.phase)
%         TB.simData2.phase = string(TB.simData2.phase);
%         TB.simData2.phase(TB.simData2.phase=="1") = "ES";
%         TB.simData2.phase(TB.simData2.phase=="2") = "EE";
%     end
end


% For each data table, overwrite OV and VD with the data from OV_2 and VD_2
TB.realData.OV = TB.realData.OV_2;
TB.realData.VD = TB.realData.VD_2;

% Rename columns in real data to match simulated data
renameMap = containers.Map({'chose_right', 'rtime', 'FirstFixDur', 'MiddleFixDur', 'FinalFixDur', 'eachMiddleFixDur'}, ...
                           {'Choice', 'rt', 'DwellFirst', 'DwellMid', 'DwellFinal', 'eachDwellMiddle'});
%chose_right
%choice
realVarNames = TB.realData.Properties.VariableNames; 

for k = 1:length(realVarNames)
    if isKey(renameMap, realVarNames{k})
        TB.realData.Properties.VariableNames{k} = renameMap(realVarNames{k});
    end
end

% variable of interst for descriptive results
VOI = {'Choice';'rt';'Nfix';'FixLocFirstCorr';'FixLocLastCorr';'DwellFirst';'DwellMid';'DwellFinal';'eachDwellMiddle'};

RTmin = 0.25;  
RTmax = 18;
TB.realData  = TB.realData(TB.realData.rt  >= RTmin & TB.realData.rt  <= RTmax, :);
TB.simData   = TB.simData( TB.simData.rt   >= RTmin & TB.simData.rt   <= RTmax, :);
% TB.simData2  = TB.simData2(TB.simData2.rt  >= RTmin & TB.simData2.rt  <= RTmax, :);


%% Descriptive statistics
switch answer
    case {'EE', 'ES'}
        % For these cases, group by OV
        for kData = 1:2
            TB_OV_Sub.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'sub_id', 'OV'});
    
            TB_Sub.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'sub_id'});
    
            TB_OV.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'OV'});
        end
    case 'ESEE'
        % For ESEE, group by phase.
        for kData = 1:3
            TB_OV_Sub.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'sub_id', 'phase'});
    
            TB_Sub.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'sub_id'});
    
            TB_OV.(DataName{kData}) = varfun(@(x)mean(x,'omitnan'), TB.(DataName{kData}), 'InputVariables', VOI, ...
                'GroupingVariables', {'phase'});
        end
end


switch answer
    case 'ESEE'                    % ES + EE pooled
        groupName   = 'phase';
        groups      = unique(TB_OV_Sub.(DataName{1}).phase);
        groupLabels = {'ES','EE'}; 
    case {'ES','EE'}               % single phase
        groupName   = 'OV';
        groups      = unique(TB_OV_Sub.(DataName{1}).OV);
        groupLabels = {'L','M','H'};   % low / medium / high OV
    otherwise
        error('Unexpected value for ANSWER: %s',answer);
end
nGroups = numel(groups);




% Figure 5 Plotting
% Parameters
fontsize = 24;
markerS = 10;
lineW   = 3;
errorCap= 10;
white   = [1 1 1];
black   = [0 0 0];
gray    = [.4 .4 .4];
colormatrix = [0 0.749 1; ...       % For real data
               0.6 0.196 0.8; ...   % for simData
               0 0.749 0];          % for simData2
colormatrix = colormatrix(realData,:);

% Variable names for plotting
variableName = {'rt (sec)', 'chose S %'};
variableIDX  = {'Fun_rt', 'Fun_Choice'};

% Create the figure
figure('Position', [10 10 1500 400]); hold on

for kData = 1:2 % Loop over datasets: 1=real, 2=simData, 3=simData2
    clear indData
    % Convert correct proportion to percentage if necessary
    TB_OV_Sub.(DataName{kData}).Fun_Choice = TB_OV_Sub.(DataName{kData}).Fun_Choice * 100;
    indData = TB_OV_Sub.(DataName{kData});
    ts = tinv(0.975, nSubj(kData)-1);  % t-score
    
    for k_var = 1:2  
        subplot(1,2,k_var); hold on
        
        if kData > 1
            tempData = cell(1, nGroups);
            for kGroup = 1:nGroups
                idx = strcmp(string(indData.(groupName)), string(groups(kGroup)));
                tempData{kGroup} = indData.(variableIDX{k_var})(idx);
            end
            inp = cellfun(@nanmean, tempData);
            err = ts .* cellfun(@nanstd, tempData) ./ sqrt(nSubj(kData));
    
            if kData == 2
                plot((1:nGroups)+0.2, inp, '-', 'LineWidth', lineW, 'MarkerSize', markerS, ...
                'MarkerEdgeColor', gray+0.3, 'Color', gray);
            else
                plot((1:nGroups)-0.2, inp, '-', 'LineWidth', lineW, 'MarkerSize', markerS, ...
                'MarkerEdgeColor', black, 'Color', black);
            end
        end
        
        for kGroup = 1:nGroups
            idx = strcmp(string(indData.(groupName)), string(groups(kGroup)));
            tempData = indData.(variableIDX{k_var})(idx);
            inp_val = nanmean(tempData);
            err_val = ts * (nanstd(tempData)) / sqrt(nSubj(kData));
            
            switch kData
                case 1  % real data
                    bar(kGroup, inp_val, 'FaceColor', colormatrix(kData,:), 'EdgeColor', white, ...
                        'LineWidth', 1.5, 'BarWidth', 0.7);
                    errorbar(kGroup, inp_val, err_val, '.-', 'LineWidth', 2, 'CapSize', 6, ...
                        'Color', colormatrix(kData,:) * 0.6);
                case 2  % simulated data (model 1)
                    plot(kGroup+0.2, inp_val, 's', 'LineWidth', 2, 'MarkerFaceColor', white, ...
                        'MarkerEdgeColor', black, 'MarkerSize', 15);
                    errorbar(kGroup+0.2, inp_val, err_val, '.', 'LineWidth', 2, 'CapSize', 8, ...
                        'Color', black);
                case 3  % simulated data (model 2)
                    plot(kGroup-0.2, inp_val, 'd', 'LineWidth', 2, 'MarkerFaceColor', gray, ...
                        'MarkerEdgeColor', black, 'MarkerSize', 15);
                    errorbar(kGroup-0.2, inp_val, err_val, '.', 'LineWidth', 2, 'CapSize', 8, ...
                        'Color', black);
            end
        end
        
        xlabel(groupName);  
        ylabel(variableName{k_var});
        xlim([0 nGroups+1]);
        xticks(1:nGroups);
        xticklabels(groupLabels);
        set(gca, 'FontSize', fontsize);
    end
end


