function combined_plotEmpiricalData(matOrCsvFile)
% combined_plotEmpiricalData loading and plotting empirical + simulated overlay

% Empirical and simulated file paths
empiricalPath = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/Data/OVParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv';
simulatedPath = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/NewFixCode_EXP2/ppc_full_bayesian/sim_1000_postdraws/OV_6_metrics_summary_mean.csv';

T = readtable(empiricalPath);
vars = T.Properties.VariableNames;

exclude = [6, 14, 20, 26, 2, 9, 18];
maskExclude = ~ismember(T.sub_id, exclude);
maskPhase = strcmp(T.phase,'ES');
T = T(maskExclude & maskPhase,:);

RTvec = T.rtime;
choiceEall = (T.cho == 1);
dwaSall = T.DwellTimeAdvantage;
rtMask = (RTvec >= 0.250) & (RTvec <= 20);

RTvec = RTvec(rtMask);
choiceEall = choiceEall(rtMask);
dwaSall = dwaSall(rtMask);
T = T(rtMask,:);

subs = unique(T.sub_id);
nS = numel(subs);
pE_byS = nan(nS,5);
pS_byS = nan(nS,5);

for i = 1:nS
    sid = subs(i);
    sel = T.sub_id == sid;
    RTs = RTvec(sel);
    chE = choiceEall(sel);
    edges = quantile(RTs,4);
    rb = discretize(RTs,[-inf edges inf]);
    for b = 1:5
        pE_byS(i,b) = mean(chE(rb==b));
        pS_byS(i,b) = mean(~chE(rb==b));
    end
end

m_qRT_E = mean(pE_byS,1); sem_qRT_E = std(pE_byS,0,1)/sqrt(nS);
m_qRT_S = mean(pS_byS,1); sem_qRT_S = std(pS_byS,0,1)/sqrt(nS);

edgesD = quantile(dwaSall,4);
db = discretize(dwaSall,[-inf edgesD inf]);
m_pS = nan(1,5); sem_pS = nan(1,5);
m_RT = nan(1,5); sem_RT = nan(1,5);
m_RT_E = nan(1,5); sem_RT_E = nan(1,5);
m_RT_S = nan(1,5); sem_RT_S = nan(1,5);

for b = 1:5
    idx = (db == b);
    N = sum(idx);
    pS = mean(~choiceEall(idx));
    m_pS(b) = pS;
    sem_pS(b) = sqrt(pS*(1 - pS)/N);
    
    theseRT = RTvec(idx);
    m_RT(b) = mean(theseRT);
    sem_RT(b) = std(theseRT)/sqrt(N);
    
    idxE = idx & choiceEall;
    idxS = idx & ~choiceEall;
    
    RT_E = RTvec(idxE);
    RT_S = RTvec(idxS);
    m_RT_E(b) = mean(RT_E);
    sem_RT_E(b) = std(RT_E)/sqrt(numel(RT_E));
    m_RT_S(b) = mean(RT_S);
    sem_RT_S(b) = std(RT_S)/sqrt(numel(RT_S));
end

m_qRT = m_qRT_E;
sem_qRT = sem_qRT_E;

% P(correct) binning for empirical
m_corrProb = NaN(1,5); sem_corrProb = NaN(1,5);
if ismember('DwelltimeAdvantageCorrect',vars) && ismember('corr',vars)
    DAC = T.DwelltimeAdvantageCorrect;
    ACC = double(T.corr);
    valid = isfinite(DAC) & isfinite(ACC);
    DAC = DAC(valid);
    ACC = ACC(valid);
    if numel(DAC) >= 5
        q = prctile(DAC,[20 40 60 80]);
        q = adjustEdges(q);
        bins = discretize(DAC,[-inf q inf]);
        for b = 1:5
            idx = bins == b;
            N = sum(idx);
            if N > 0
                p = mean(ACC(idx));
                m_corrProb(b) = p;
                sem_corrProb(b) = sqrt(p*(1-p)/N);
            end
        end
    end
end

% loading simulated summary
simT = readtable(simulatedPath);

s_qRT = simT.Mean_qRT_E(:)';
s_qRT_sem = simT.SEM_qRT_E(:)';
s_qRT_S = simT.Mean_qRT_S(:)';
s_qRT_S_sem = simT.SEM_qRT_S(:)';
s_pS = simT.Mean_pS(:)'; s_pS_sem = simT.SEM_pS(:)';
s_RT = simT.Mean_RT(:)'; s_RT_sem = simT.SEM_RT(:)';
s_corrProb = simT.Mean_CorrectProb(:)'; s_corrProb_sem = simT.SEM_CorrectProb(:)';

% PLotting
combined_simulationESmodel_plotWithSEM(...
    m_qRT, sem_qRT, ...
    m_qRT_S, sem_qRT_S, ...
    m_pS, sem_pS, ...
    m_RT, sem_RT, ...
    m_corrProb, sem_corrProb, ...
    s_qRT, s_qRT_sem, ...
    s_qRT_S, s_qRT_S_sem, ...
    s_pS, s_pS_sem, ...
    s_RT, s_RT_sem, ...
    s_corrProb, s_corrProb_sem);

sgtitle('Empirical (bars) + Simulated (dotted)');

end

function q = adjustEdges(q)
  epsStep = max(1, range(q)) * 1e-12;
  for k = 2:numel(q)
      if q(k) <= q(k-1), q(k) = q(k-1) + epsStep; end
  end
end
