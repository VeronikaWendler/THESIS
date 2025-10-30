function plotEmpiricalData(matOrCsvFile)
% plotEmpiricalData  Load CSV or MAT, plot, save 
% There are differences in how some columns were saved:
% CSV:  sub_id, phase, rtime, cho(1=E,2=S), DwellTimeAdvantage (this is S-E )
% SIM:  sub_id, (no phase), RT, Choice(0=E,1=S), DwellDiff (S-E same as R-L)

  % choose file
  if nargin<1 || isempty(matOrCsvFile)
    dataDir = fullfile(pwd,'Sim');
    files   = [dir(fullfile(dataDir,'*.csv')); dir(fullfile(dataDir,'*.mat'))];
    if isempty(files), error('No CSV or MAT found in %s',dataDir); end
    if numel(files)>1
      [f,p] = uigetfile({'*.csv;*.mat','Data files (*.csv,*.mat)'}, dataDir);
      if f==0, return; end
      matOrCsvFile = fullfile(p,f);
    else
      matOrCsvFile = fullfile(dataDir,files(1).name);
    end
  end

  [~,~,ext] = fileparts(matOrCsvFile);
  switch lower(ext)
    case '.csv'
      T = readtable(matOrCsvFile);
    case '.mat'
      T = loadTableFromMat(matOrCsvFile);
    otherwise
      error('Unsupported extension %s',ext);
  end
  
  % aggregated summary CSV (means + SEMs)
  vars = T.Properties.VariableNames;
  isSummary = all(ismember({ ...
      'Bin', ...
      'Mean_qRT_E','SEM_qRT_E', ...
      'Mean_qRT_S','SEM_qRT_S', ...
      'Mean_pS','SEM_pS', ...
      'Mean_RT','SEM_RT', ...
      'Mean_RT_E','SEM_RT_E', ...
      'Mean_RT_S','SEM_RT_S', ...
      'Mean_CorrectProb','SEM_CorrectProb'}, vars));

  if isSummary
      m_qRT_E   = T.Mean_qRT_E(:)';   sem_qRT_E   = T.SEM_qRT_E(:)';
      m_qRT_S   = T.Mean_qRT_S(:)';   sem_qRT_S   = T.SEM_qRT_S(:)';
      m_pS      = T.Mean_pS(:)';      sem_pS      = T.SEM_pS(:)';
      m_RT      = T.Mean_RT(:)';      sem_RT      = T.SEM_RT(:)';
      m_RT_E    = T.Mean_RT_E(:)';    sem_RT_E    = T.SEM_RT_E(:)';
      m_RT_S    = T.Mean_RT_S(:)';    sem_RT_S    = T.SEM_RT_S(:)';
      m_corrProb= T.Mean_CorrectProb(:)'; sem_corrProb = T.SEM_CorrectProb(:)';
      m_qRT   = m_qRT_E;
      s_qRT   = sem_qRT_E;

      simulationESmodel_plotWithSEM(NaN, NaN, NaN, ...
          m_qRT,   s_qRT, ...
          m_pS,    sem_pS, ...
          m_RT,    sem_RT, ...
          m_corrProb, sem_corrProb);

      [~,baseName,~] = fileparts(matOrCsvFile);
      sgtitle(sprintf('Summary means from %s', baseName), 'Interpreter','none');

      save([baseName '_plotted.mat'], ...
          'm_qRT_E','sem_qRT_E','m_qRT_S','sem_qRT_S', ...
          'm_pS','sem_pS','m_RT','sem_RT', ...
          'm_RT_E','sem_RT_E','m_RT_S','sem_RT_S', ...
          'm_corrProb','sem_corrProb');

      return
  end



  % figure out format by variable names
  vars = T.Properties.VariableNames;
  isReal = ismember('cho',vars);      % in my empirical (real) CSV data this is the name
  isSim  = ismember('Choice',vars);   % simulated MAT (similar to Chih-Chung) 


  % filter for subjects and ES phase 
  exclude = [1,4,5,6,14,99];       % for EXP2 ..6, 14, 20, 26, 2, 9, 18  for EXP1: 1,4,5,6,14,99
  maskExclude = ~ismember(T.sub_id,exclude);

  if ismember('phase',vars)
      maskPhase = strcmp(T.phase,'ES');
  else
      maskPhase = true(size(maskExclude)); % sim files already ES-only
  end

% if interested in specific fixation constraints 
%   maskFix = true(height(T),1);
%   if strcmpi(ext,'.csv')
%       needCols = {'RightFixNR','LeftFixNR'};
%       if all(ismember(needCols, vars))
%           maskFix = (T.RightFixNR > 0) & (T.LeftFixNR > 0);
%       else
%           warning('skipping fix-count filter');
%       end
%   end
% 
%   T = T(maskExclude & maskPhase & maskFix, :);

  T = T(maskExclude & maskPhase,:);

  if isReal
      RTvec      = T.rtime;
      choiceEall = (T.cho==1);                  % 1 = E chosen (in my empirical csv data, cho is 2 for S and 1 for E chosen)
      dwaSall    = T.DwellTimeAdvantage;        % S - E
  else % sim
      RTvec      = T.rt;
      choiceEall = (T.Choice==0);               % do this: Choice == 1 if cdoe involves E as uppper bound, otherwise set to 0 
      dwaSall    = T.DwellDiff;                 % S - E
  end
 
  %RT restriction criteria
  rtMask = (RTvec >= 0) & (RTvec <= 18);

  RTvec      = RTvec(rtMask);
  choiceEall = choiceEall(rtMask);
  dwaSall    = dwaSall(rtMask);
  T          = T(rtMask,:);

  % analysis - unchanged (this is similar to Sebastian's code)
  subs     = unique(T.sub_id);
  nS       = numel(subs);
  pE_byS   = nan(nS,5);
  pS_byS   = nan(nS,5);

  for i = 1:nS
    sid  = subs(i);
    sel  = T.sub_id==sid;
    RTs  = RTvec(sel);
    chE  = choiceEall(sel);
    edges = quantile(RTs,4);
    rb    = discretize(RTs,[-inf edges inf]);
    for b=1:5
      pE_byS(i,b) = mean(   chE(rb==b) );
      pS_byS(i,b) = mean( ~ chE(rb==b) );
    end
  end

  m_qRT_E   = mean(pE_byS,1);    sem_qRT_E = std(pE_byS,0,1)/sqrt(nS);
  m_qRT_S   = mean(pS_byS,1);    sem_qRT_S = std(pS_byS,0,1)/sqrt(nS);

  edgesD = quantile(dwaSall,4);
  db     = discretize(dwaSall,[-inf edgesD inf]);
  m_pS   = nan(1,5); sem_pS = nan(1,5);
  m_RT   = nan(1,5); sem_RT = nan(1,5);
  m_RT_E = nan(1,5); sem_RT_E = nan(1,5);
  m_RT_S = nan(1,5); sem_RT_S = nan(1,5);

  for b=1:5
    idx        = (db==b);
    N          = sum(idx);
    pS         = mean(~choiceEall(idx));
    m_pS(b)    = pS;
    sem_pS(b)  = sqrt(pS*(1-pS)/N);

    theseRT    = RTvec(idx);
    m_RT(b)    = mean(theseRT);
    sem_RT(b)  = std(theseRT)/sqrt(N);

    idxE       = idx &  choiceEall;
    idxS       = idx & ~choiceEall;

    RT_E       = RTvec(idxE);
    m_RT_E(b)  = mean(RT_E);
    sem_RT_E(b)= std(RT_E)/sqrt(numel(RT_E));

    RT_S       = RTvec(idxS);
    m_RT_S(b)  = mean(RT_S);
    sem_RT_S(b)= std(RT_S)/sqrt(numel(RT_S));
  end

% for simulation and csv
m_qRT   = m_qRT_E;
sem_qRT = sem_qRT_E;
% for plot 4 P(correct) by dwell advantage for correct option
m_corrProb   = NaN(1,5);
sem_corrProb = NaN(1,5);

if isReal
    DAC = T.DwelltimeAdvantageCorrect;   % >0 = more dwell on correct
    ACC = double(T.corr);                % 0/1 accuracy

    valid = isfinite(DAC) & isfinite(ACC);
    DAC   = DAC(valid);
    ACC   = ACC(valid);

    if numel(DAC) >= 5
        q = prctile(DAC,[20 40 60 80]);

        % avoid duplicate edges (flat distributions)
        epsStep = max(1, range(DAC)) * 1e-12;
        for k = 2:numel(q)
            if q(k) <= q(k-1), q(k) = q(k-1) + epsStep; end
        end

        edges = [-inf q inf];
        bins  = discretize(DAC, edges);

        for b = 1:5
            idx = (bins == b);
            N   = sum(idx);
            if N > 0
                p  = mean(ACC(idx));
                m_corrProb(b)   = p;
                sem_corrProb(b) = sqrt(p*(1-p)/N);
            end
        end
    end

else
    % SIM path
    if all(ismember({'DwellDiff','Correct','Vl','Vr'}, vars))
        valid = ~isnan(T.DwellDiff) & (T.Vl ~= T.Vr);
        DAC   = sign(T.Vr(valid) - T.Vl(valid)) .* T.DwellDiff(valid);
        ACC   = T.Correct(valid);

        if numel(DAC) >= 5
            q = prctile(DAC,[20 40 60 80]);
            epsStep = max(1, range(DAC))*1e-12;
            for k = 2:numel(q)
                if q(k) <= q(k-1), q(k) = q(k-1) + epsStep; end
            end
            edges = [-inf q inf];
            bins  = discretize(DAC, edges);

            for b = 1:5
                idx = (bins == b);
                N   = sum(idx);
                if N > 0
                    p  = mean(ACC(idx));
                    m_corrProb(b)   = p;
                    sem_corrProb(b) = sqrt(p*(1-p)/N);
                end
            end
        end
    end
end


  % callign  the plotting function
  simulationESmodel_plotWithSEM(NaN, NaN, NaN, ...
    m_qRT, sem_qRT, ...
    m_pS, sem_pS, ...
    m_RT, sem_RT, ...
    m_corrProb, sem_corrProb);


  sgtitle('Empirical (phase=ES, excl subs [1,4,5,6,14,99])');

  save('empirical_metrics.mat', ...
    'm_qRT_E','sem_qRT_E', ...
    'm_qRT_S','sem_qRT_S', ...
    'm_pS','sem_pS', ...
    'm_RT','sem_RT', ...
    'm_RT_E','sem_RT_E', ...
    'm_RT_S','sem_RT_S');


% table for CSV export
summaryStats = table( ...
    (1:5)', ...
    m_qRT_E', sem_qRT_E', ...
    m_qRT_S', sem_qRT_S', ...
    m_pS', sem_pS', ...
    m_RT', sem_RT', ...
    m_RT_E', sem_RT_E', ...
    m_RT_S', sem_RT_S', ...
    m_corrProb', sem_corrProb', ...
    'VariableNames', { ...
        'Bin', ...
        'Mean_qRT_E', 'SEM_qRT_E', ...
        'Mean_qRT_S', 'SEM_qRT_S', ...
        'Mean_pS', 'SEM_pS', ...
        'Mean_RT', 'SEM_RT', ...
        'Mean_RT_E', 'SEM_RT_E', ...
        'Mean_RT_S', 'SEM_RT_S', ...
        'Mean_CorrectProb', 'SEM_CorrectProb' ...
    });

writetable(summaryStats, 'For_paper_7_metrics_summary_sim10.csv');


end


% Some functions
%%----------------------------------------------------%%
function T = loadTableFromMat(fname)
  S = load(fname);
  if isfield(S,'TB') && isfield(S.TB,'Garcia') && istable(S.TB.Garcia)
      T = S.TB.Garcia; return;
  end
  if isfield(S,'TBsim') && istable(S.TBsim)
      T = S.TBsim; return;
  end
  fns = fieldnames(S);
  for k = 1:numel(fns)
      T = findFirstTable(S.(fns{k}));
      if ~isempty(T), return; end
  end
  error('No table in %s', fname);
end

function T = findFirstTable(x)
  T = [];
  if istable(x), T = x; return; end
  if isstruct(x)
      f = fieldnames(x);
      for i=1:numel(f)
          T = findFirstTable(x.(f{i}));
          if ~isempty(T), return; end
      end
  elseif iscell(x)
      for i=1:numel(x)
          T = findFirstTable(x{i});
          if ~isempty(T), return; end
      end
  end
end


%%----------------------------------------------------%%


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for simulation data %%%%%%%%%%%%%%%%%%%
%   % rename for plotting
%   m_qRT   = m_qRT_E;
%   sem_qRT = sem_qRT_E;

%   % Accuracy vs Dwell Time Advantage for Correct Option (assumes Vl, Vr) --
%   m_corrProb  = NaN(1,5);
%   sem_corrProb= NaN(1,5);
% 
%   if isSim && ismember('DwellDiff', T.Properties.VariableNames) && ismember('Correct', T.Properties.VariableNames)
%     
%       if ~all(ismember({'Vl','Vr'}, T.Properties.VariableNames))
%           warning('Skipping plot %s', strjoin(T.Properties.VariableNames, ', '));
%       else
%      
%           valid = ~isnan(T.DwellDiff) & (T.Vl ~= T.Vr);
% 
%           % dwell-time advantage for the correct option
%           DAC = sign(T.Vr(valid) - T.Vl(valid)) .* T.DwellDiff(valid);  % >0  dwell favors higher-value option
%           ACC = T.Correct(valid);                                       % 1 = chose higher value
% 
%           % need enough data
%           if numel(DAC) >= 5
%               % 5 bins quantiles
%               q = prctile(DAC,[20 40 60 80]);
%               epsStep = max(1, range(DAC))*1e-12;
%               for k = 2:numel(q)
%                   if q(k) <= q(k-1), q(k) = q(k-1) + epsStep; end
%               end
%               edges = [-inf q inf];
%               bins  = discretize(DAC, edges);
% 
%               for b = 1:5
%                   idx = (bins == b);
%                   N   = sum(idx);
%                   if N > 0
%                       p  = mean(ACC(idx));
%                       m_corrProb(b)   = p;
%                       sem_corrProb(b) = sqrt(p*(1-p)/N);
%                   end
%               end
%           else
%               warning('Not enough trials for plot n=%d).', numel(DAC));
%           end
%       end
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% 
% %% harmonize columns
%   if isReal
%       RTvec      = T.rtime;
%       choiceEall = (T.cho==1);                  % 1 = E chosen
%       if ismember('DwellTimeAdvantage',vars)
%           dwaSall = T.DwellTimeAdvantage;       % S - E
%       else
%           error('Real data missing DwellTimeAdvantage column.');
%       end
%   else % sim
%       % be robust to 'rt' or 'RT'
%       if ismember('rt',vars), RTvec = T.rt; elseif ismember('RT',vars), RTvec = T.RT; else, error('Sim data missing rt/RT.'); end
%       choiceEall = (T.Choice==0);               % 0 = E chosen
%       if ismember('DwellDiff',vars)
%           dwaSall = T.DwellDiff;                % S - E
%       else
%           error('Sim data missing DwellDiff column.');
%       end
%   end
% 
%   % clip RTs to [0,8]
%   rtMask = (RTvec >= 0) & (RTvec <= 8);
%   RTvec      = RTvec(rtMask);
%   choiceEall = choiceEall(rtMask);
%   dwaSall    = dwaSall(rtMask);
%   T          = T(rtMask,:);
% 
%   %% ----- ANALYSIS -----
%   nBins   = 5;
%   qProbs  = [0.2 0.4 0.6 0.8]; % explicit quantiles for 5 bins
%   subs    = unique(T.sub_id);
%   nS      = numel(subs);
% 
%   % ---- RT-QUANTILE ANALYSIS (within-subject first, as before) ----
%   pE_byS = nan(nS,nBins);
%   pS_byS = nan(nS,nBins);
% 
%   for i = 1:nS
%     sid   = subs(i);
%     sel   = (T.sub_id==sid);
%     RTs   = RTvec(sel);
%     chE   = choiceEall(sel);
% 
%     % subject-specific RT quantiles (5 bins via 20/40/60/80 cutpoints)
%     edgesRT = quantile(RTs, qProbs);
%     rb      = discretize(RTs, [-inf edgesRT inf]);
% 
%     % per-bin choice proportions for this subject
%     for b = 1:nBins
%       pE_byS(i,b) = mean( chE(rb==b), 'omitnan' );
%       pS_byS(i,b) = mean(~chE(rb==b), 'omitnan' );
%     end
%   end
% 
%   m_qRT_E   = mean(pE_byS,1,'omitnan');
%   sem_qRT_E = std(pE_byS,0,1,'omitnan')/sqrt(nS);
% 
%   m_qRT_S   = mean(pS_byS,1,'omitnan');
%   sem_qRT_S = std(pS_byS,0,1,'omitnan')/sqrt(nS);
% 
%   % for plotting (legacy function expects this)
%   m_qRT     = m_qRT_E;
%   sem_qRT   = sem_qRT_E;
% 
%   % ---- DWELL-TIME-ADVANTAGE ANALYSIS (NOW within-subject first) ----
%   % We compute, for each subject and each dwell-time bin:
%   %   - pS (prob choose S)
%   %   - mean RT (all trials)
%   %   - mean RT_E (trials where E chosen)
%   %   - mean RT_S (trials where S chosen)
%   pS_dwell_byS   = nan(nS,nBins);
%   RT_dwell_byS   = nan(nS,nBins);
%   RT_E_dwell_byS = nan(nS,nBins);
%   RT_S_dwell_byS = nan(nS,nBins);
% 
%   for i = 1:nS
%     sid  = subs(i);
%     sel  = (T.sub_id==sid);
% 
%     dwa  = dwaSall(sel);   % S−E dwell advantage for this subject
%     rts  = RTvec(sel);
%     chE  = choiceEall(sel);
% 
%     % subject-specific dwell quantiles (5 bins via 20/40/60/80 cutpoints)
%     edgesD = quantile(dwa, qProbs);
%     db     = discretize(dwa, [-inf edgesD inf]);
% 
%     for b = 1:nBins
%       idx   = (db==b);
%       pS_dwell_byS(i,b)   = mean(~chE(idx), 'omitnan');   % prob choose S
%       RT_dwell_byS(i,b)   = mean(  rts(idx), 'omitnan');
% 
%       idxE  = idx &  chE;     % E chosen
%       idxS  = idx & ~chE;     % S chosen
%       RT_E_dwell_byS(i,b) = mean(rts(idxE), 'omitnan');
%       RT_S_dwell_byS(i,b) = mean(rts(idxS), 'omitnan');
%     end
%   end
% 
%   % group means and SEM across subjects
%   m_pS     = mean(pS_dwell_byS,1,'omitnan');
%   sem_pS   = std( pS_dwell_byS,0,1,'omitnan')/sqrt(nS);
% 
%   m_RT     = mean(RT_dwell_byS,1,'omitnan');
%   sem_RT   = std( RT_dwell_byS,0,1,'omitnan')/sqrt(nS);
% 
%   m_RT_E     = mean(RT_E_dwell_byS,1,'omitnan');
%   sem_RT_E   = std( RT_E_dwell_byS,0,1,'omitnan')/sqrt(nS);
% 
%   m_RT_S     = mean(RT_S_dwell_byS,1,'omitnan');
%   sem_RT_S   = std( RT_S_dwell_byS,0,1,'omitnan')/sqrt(nS);
% 
%   % ---- PLOT (uses your existing plotting helper) ----
%   simulationESmodel_plotWithSEM(NaN,NaN,NaN, ...
%       m_qRT,sem_qRT, m_pS,sem_pS, m_RT,sem_RT);
%   sgtitle('Empirical (phase=ES, excl subs [1,4,5,6,14,99]) – dwell is within-subject');
% 
%   % ---- SAVE METRICS ----
%   save('empirical_metrics.mat', ...
%     'm_qRT_E','sem_qRT_E', ...
%     'm_qRT_S','sem_qRT_S', ...
%     'm_pS','sem_pS', ...
%     'm_RT','sem_RT', ...
%     'm_RT_E','sem_RT_E', ...
%     'm_RT_S','sem_RT_S');
% end
% 
% %%----------------------------------------------------%%
% function T = loadTableFromMat(fname)
% % Try the known fields, then recurse
%   S = load(fname);
%   if isfield(S,'TB') && isfield(S.TB,'Garcia') && istable(S.TB.Garcia)
%       T = S.TB.Garcia; return;
%   end
%   if isfield(S,'TBsim') && istable(S.TBsim)
%       T = S.TBsim; return;
%   end
%   % fallback recursive search
%   fns = fieldnames(S);
%   for k = 1:numel(fns)
%       T = findFirstTable(S.(fns{k}));
%       if ~isempty(T), return; end
%   end
%   error('No recognizable table in %s', fname);
% end
% 
% function T = findFirstTable(x)
%   T = [];
%   if istable(x), T = x; return; end
%   if isstruct(x)
%       f = fieldnames(x);
%       for i=1:numel(f)
%           T = findFirstTable(x.(f{i}));
%           if ~isempty(T), return; end
%       end
%   elseif iscell(x)
%       for i=1:numel(x)
%           T = findFirstTable(x{i});
%           if ~isempty(T), return; end
%       end
%   end
% end
% 




% function plotEmpiricalData()
% % plotEmpiricalData  Load CSV, draw the 3‐panel plot from @Gluth, save
% % empirical_metrics.mat with data 
% % needed in csv: sub_id, phase (LE;ES;EE;SP), rtime, cho, DwellTimeAdvantage
% 
%   dataDir = fullfile(pwd,'data_BGarcia');
%   csvs    = dir(fullfile(dataDir,'*.csv'));
%   if isempty(csvs)
%     error('Error!', dataDir);
%   elseif numel(csvs)>1
%     [f,p] = uigetfile('*.csv', dataDir);
%     if f==0, return; end
%     csvFile = fullfile(p,f);
%   else
%     csvFile = fullfile(dataDir,csvs(1).name);
%   end
% 
%   T = readtable(csvFile);
%   exclude = [1,4,5,6,14,99];                         % just for EXP1 - 1,4,5,6,14,99  for EXP2 - 6,14,20,26,2,9,18
%   T = T(~ismember(T.sub_id,exclude) & strcmp(T.phase,'ES'),:);
% 
%   RTvec      = T.rtime;                % trial RT
%   choiceEall = (T.cho==1);             % true=E chosen
%   dwaSall    = T.DwellTimeAdvantage;   % dwell S–E
% 
%   % P(E) by RT‐quantile subject‐level bins 
%   subs     = unique(T.sub_id);
%   nS       = numel(subs);
%   pE_byS   = nan(nS,5);
%   pS_byS   = nan(nS,5);     
% 
%   for i = 1:nS
%       sid  = subs(i);
%       sel  = T.sub_id==sid;
%       RTs  = RTvec(sel);
%       chE  = choiceEall(sel);
%       edges = quantile(RTs,4);  
%       rb    = discretize(RTs,[-inf edges inf]);
%       for b = 1:5
%           pE_byS(i,b) = mean(   chE(rb==b)    );   % P(E)
%           pS_byS(i,b) = mean(~ chE(rb==b)    );   % P(S) = 1 − P(E)
%       end
%   end
% 
%   % mean across subs
%   m_qRT_E   = mean(pE_byS,1);
%   sem_qRT_E = std (pE_byS,0,1)/sqrt(nS);
%   m_qRT_S   = mean(pS_byS,1);
%   sem_qRT_S = std (pS_byS,0,1)/sqrt(nS);
% 
%   % DWA bins 
%   edgesD = quantile(dwaSall,4);
%   db     = discretize(dwaSall,[-inf edgesD inf]);
%   m_pS   = nan(1,5); sem_pS = nan(1,5);
%   m_RT   = nan(1,5); sem_RT = nan(1,5);
%   for b=1:5
%     idx = (db==b);
%     N   = sum(idx);
%     % P(choose S)
%     pS       = mean(~choiceEall(idx));
%     m_pS(b)  = pS;
%     sem_pS(b)= sqrt(pS*(1-pS)/N);
%     % RT
%     theseRT  = RTvec(idx);
%     m_RT(b)  = mean(theseRT);
%     sem_RT(b)= std(theseRT)/sqrt(N);
% 
% 
%     idxE       = idx &  choiceEall;  % trials where E was chosen
%     idxS       = idx & ~choiceEall;  % trials where S was chosen
% 
%     RT_E       = RTvec(idxE);
%     m_RT_E(b)  = mean(RT_E);
%     sem_RT_E(b)= std(RT_E)/sqrt(numel(RT_E));
% 
%     RT_S       = RTvec(idxS);
%     m_RT_S(b)  = mean(RT_S);
%     sem_RT_S(b)= std(RT_S)/sqrt(numel(RT_S));
% 
% 
%   end
% 
%   % renaming for plotting
%   m_qRT   = m_qRT_E;
%   sem_qRT = sem_qRT_E;
% 
%   % plot function
%   simulationESmodel_plotWithSEM(NaN,NaN,NaN, ...
%     m_qRT,sem_qRT, m_pS,sem_pS, m_RT,sem_RT);
%   sgtitle('Empirical (phase=ES, excl subs [1,4,5,6,14,99])');
% 
%  save('empirical_metrics.mat', ...
%   'm_qRT_E','sem_qRT_E', ...   
%   'm_qRT_S','sem_qRT_S', ...   
%   'm_pS','sem_pS', ...         
%   'm_RT','sem_RT',...
%   'm_RT_E','sem_RT_E',...
%   'm_RT_S','sem_RT_S');
% end




% function plotEmpiricalData()
%   % 1) load & filter
%   dataDir = fullfile(pwd,'data');
%   csvs    = dir(fullfile(dataDir,'*.csv'));
%   if numel(csvs)>1
%     [f,p] = uigetfile('*.csv','Select CSV',dataDir); 
%     if f==0, return; end
%     csvFile = fullfile(p,f);
%   else
%     csvFile = fullfile(dataDir,csvs(1).name);
%   end
%   T = readtable(csvFile);
%   exclude = [1,4,5,6,14,99];
%   T = T(~ismember(T.sub_id,exclude) & strcmp(T.phase,'ES'),:);
% 
%   RTvec      = T.rtime;
%   choiceEall = (T.cho==1);
%   dwaSall    = T.DwellTimeAdvantage;
% 
%   %% 2) RT‐quantile panel: **subject‐level** bins
%   subs = unique(T.sub_id);
%   nS   = numel(subs);
%   pE_byS = nan(nS,5);
%   for i = 1:nS
%     sid      = subs(i);
%     sel      = T.sub_id==sid;
%     RTsub    = RTvec(sel);
%     choiceES = choiceEall(sel);
% 
%     edgesRT  = quantile(RTsub,4);
%     rtBin    = discretize(RTsub,[-inf edgesRT inf]);
% 
%     for b = 1:5
%       pE_byS(i,b) = mean( choiceES(rtBin==b) );
%     end
%   end
% 
%   m_qRT   = mean(pE_byS,1);
%   sem_qRT = std(pE_byS,0,1)/sqrt(nS);
% 
%   %% 3) DWA panels: **global** bins (just like your sims)
%   edgesD  = quantile(dwaSall,4);
%   dwaBin  = discretize(dwaSall,[-inf edgesD inf]);
% 
%   m_pS   = nan(1,5); sem_pS = nan(1,5);
%   m_RT   = nan(1,5); sem_RT = nan(1,5);
%   for b = 1:5
%     idx = dwaBin==b;
%     N   = sum(idx);
% 
%     % P(choose S)
%     pS = mean(~choiceEall(idx));
%     m_pS(b)   = pS;
%     sem_pS(b) = sqrt(pS*(1-pS)/N);
% 
%     % RT
%     theseRT   = RTvec(idx);
%     m_RT(b)   = mean(theseRT);
%     sem_RT(b) = std(theseRT)/sqrt(N);
%   end
% 
%   %% 4) hand off to your plot helper
%   simulationESmodel_plotWithSEM( NaN,NaN,NaN, ...
%      m_qRT,sem_qRT, m_pS,sem_pS, m_RT,sem_RT );
% 
%   sgtitle(sprintf('Empirical (ES phase, excl subs %s)', mat2str(exclude)));
% end
