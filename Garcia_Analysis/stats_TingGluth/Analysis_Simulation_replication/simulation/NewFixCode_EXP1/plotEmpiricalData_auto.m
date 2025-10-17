function plotEmpiricalData_auto(matOrCsvFile, autoMode)
% plotEmpiricalData  Load CSV or MAT, plot, save 
%
% Usage:
%   plotEmpiricalData_auto()          % manual selection
%   plotEmpiricalData_auto('.mat')    % single file
%   plotEmpiricalData_auto('sim_10_postdraws2', true)   % all files in folder

  if nargin < 2
      autoMode = false;
  end

  if autoMode && isfolder(matOrCsvFile)
      dataDir = matOrCsvFile;
      files   = dir(fullfile(dataDir,'*.mat'));
      if isempty(files)
          error('No .mat files found in %s', dataDir); 
      end

      for iFile = 1:numel(files)
          fpath = fullfile(dataDir, files(iFile).name);
          fprintf('Processing %s (%d/%d)\n', files(iFile).name, iFile, numel(files));

          % recursive call for single-file mode
          plotEmpiricalData_auto(fpath, false);
      end
      return
  end

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

  % laod data
  [folder,baseName,ext] = fileparts(matOrCsvFile);
  switch lower(ext)
    case '.csv'
      T = readtable(matOrCsvFile);
    case '.mat'
      T = loadTableFromMat(matOrCsvFile);
    otherwise
      error('Unsupported extension %s',ext);
  end

  % Analysis
  vars = T.Properties.VariableNames;
  isReal = ismember('cho',vars);   % empirical CSV
  isSim  = ismember('Choice',vars);% simulation MAT

  exclude = [1,4,5,6,14,99];
  maskExclude = ~ismember(T.sub_id,exclude);
  if ismember('phase',vars)
      maskPhase = strcmp(T.phase,'ES');
  else
      maskPhase = true(size(maskExclude));
  end
  T = T(maskExclude & maskPhase,:);

  if isReal
      RTvec      = T.rtime;
      choiceEall = (T.cho==1);                  % 1 in emp data
      dwaSall    = T.DwellTimeAdvantage;        
  else % sim
      RTvec      = T.rt;
      choiceEall = (T.Choice==0);                % do this: Choice == 1 if cdoe involves E as upper bound, otherwise set to 0 
      dwaSall    = T.DwellDiff;                  % S-E    
  end

  rtMask = (RTvec >= 0) & (RTvec <= 18);
  RTvec      = RTvec(rtMask);
  choiceEall = choiceEall(rtMask);
  dwaSall    = dwaSall(rtMask);
  T          = T(rtMask,:);

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

  m_qRT   = m_qRT_E;
  sem_qRT = sem_qRT_E;
  m_corrProb   = NaN(1,5);
  sem_corrProb = NaN(1,5);

  if isReal
      DAC = T.DwelltimeAdvantageCorrect;
      ACC = double(T.corr);
      valid = isfinite(DAC) & isfinite(ACC);
      DAC   = DAC(valid);
      ACC   = ACC(valid);
      if numel(DAC) >= 5
          q = prctile(DAC,[20 40 60 80]);
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

  % Plot 
  simulationESmodel_plotWithSEM(NaN, NaN, NaN, ...
    m_qRT, sem_qRT, ...
    m_pS, sem_pS, ...
    m_RT, sem_RT, ...
    m_corrProb, sem_corrProb);

  sgtitle(sprintf('Results from %s', baseName),'Interpreter','none');

  % Save participant-level outputs 
  subjData = [];
  for i = 1:nS
      for b = 1:5
          subjData = [subjData; {subs(i), b, pE_byS(i,b), pS_byS(i,b)}];
      end
  end

  subjTable = cell2table(subjData, ...
      'VariableNames', {'Subject','Bin','ProbE','ProbS'});

  % save both summary and participant-level
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

  summaryFile = fullfile(folder, [baseName '_summary.csv']);
  subjFile    = fullfile(folder, [baseName '_subjectLevel.csv']);
  writetable(summaryStats, summaryFile);
  writetable(subjTable, subjFile);
  savefig(gcf, fullfile(folder, [baseName '_plot.fig']));
  saveas(gcf, fullfile(folder, [baseName '_plot.png']));

end

% Helpers
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
