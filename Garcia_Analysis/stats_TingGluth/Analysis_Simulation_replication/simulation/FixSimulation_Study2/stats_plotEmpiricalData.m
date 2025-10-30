function stats_plotEmpiricalData()

empiricalPath = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/Data/OVParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv';
simDir        = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/ppc_full_bayesian/sim_1000_postdraws2';
nDraws        = 1000;

excludeSubs = [6, 14, 20, 26, 2, 9, 18];

T = readtable(empiricalPath);

maskExclude = ~ismember(T.sub_id, excludeSubs);
maskPhase   = strcmp(T.phase,'ES');
rtMask      = (T.rtime >= 0) & (T.rtime <= 18);
T = T(maskExclude & maskPhase & rtMask,:);

subs = unique(T.sub_id);
nS   = numel(subs);

% empirical choice probability by RT quintile
pE_emp = nan(nS,5);
pS_emp = nan(nS,5);
for i = 1:nS
    sid = subs(i);
    sel = T.sub_id == sid;
    RTs = T.rtime(sel);
    chE = (T.cho(sel) == 1); % 1=E chosen in empirical CSV

    edges = quantile(RTs,4);
    rb = discretize(RTs,[-inf edges inf]);
    for b = 1:5
        pE_emp(i,b) = mean(chE(rb==b));
        pS_emp(i,b) = mean(~chE(rb==b));
    end
end

m_empE  = mean(pE_emp,1); sem_empE = std(pE_emp,0,1)/sqrt(nS);

pE_sim = nan(nDraws,nS,5);
pS_sim = nan(nDraws,nS,5);

for d = 1:nDraws
    fname = fullfile(simDir, sprintf('Sim_HDDM_m6_draw%03d.mat', d));
    S = load(fname); 
    TB = S.TBsim; 

    for i = 1:nS
        sid = subs(i);
        sel = TB.sub_id == sid;

        if ismember('rtime', TB.Properties.VariableNames)
            RTs = TB.rtime(sel);
            chE = (TB.cho(sel) == 1); % empirical 
        elseif ismember('rt', TB.Properties.VariableNames)
            RTs = TB.rt(sel);
            chE = (TB.Choice(sel) == 0); % sim: 0 = E, 1 = S
        else
            error('No RT column found in dataset %s', fname);
        end

        edges = quantile(RTs,4);
        rb = discretize(RTs,[-inf edges inf]);
        for b = 1:5
            pE_sim(d,i,b) = mean(chE(rb==b));
            pS_sim(d,i,b) = mean(~chE(rb==b));
        end
    end
end

% long from table for linear mixed model in R 
empRows = [];
for i = 1:nS
    for b = 1:5
        empRows = [empRows; {subs(i), b, pE_emp(i,b), 'empirical', NaN}];
    end
end

simRows = [];
for d = 1:nDraws
    for i = 1:nS
        for b = 1:5
            simRows = [simRows; {subs(i), b, pE_sim(d,i,b), 'simulated', d}];
        end
    end
end

allData = cell2table([empRows; simRows], ...
    'VariableNames', {'subject','bin','probE','source','draw'});

writetable(allData,'allData_forLMM.csv');

% PPC
meanSimE = squeeze(mean(pE_sim,2)); 
ppcMeanE = mean(meanSimE,1);
ppcCIE   = prctile(meanSimE,[2.5 97.5],1);

% Plot
figure; hold on;
errorbar(1:5, m_empE, sem_empE, 'ko-', 'LineWidth',1.5, 'MarkerFaceColor','k');
for b = 1:5
    plot([b b], ppcCIE(:,b), 'b-', 'LineWidth',4);
end
plot(1:5, ppcMeanE, 'bo--','LineWidth',1.5);
xlabel('RT Quintile'); ylabel('P(E choice)');
legend('Empirical mean ± SEM','Simulated 95% CI','Simulated mean','Location','best');
title('Posterior Predictive Check: Choice probability by RT quintile');
grid on;

end
