% Modification of Script from Chi-Chung Ting
% - Main script used to perform evidence accumulation for each trial
% function needs to be manually adjusted depending on the model

function [Choice, RT, E, tempEyeData, FixAaLL, sumdvALL, FixDur]= EvidenceAccumulate_S_upper(OV,Vl,Vr,params)

a     = params(1).*0.5;        % boundary seperation 
ndt   = params(2);             % non-decision time
z     = params(3);
beta2 = params(4)./1000;       % attention
beta3 = params(5)./1000;       % inattention E
%beta4 = params(6)./1000;       % inattention S

noisy = 0.025;                  %params(6);       % do 0.03
%theta = beta3./beta2;

% z 
spbias = ((2*z - 1)*a);       % if no bias, then z=0.5, and spbias simply 0


%% pre-determine fixations
pre_Nfix = 4; % can also assume 70 fixations before decision
FixR     = rand()> 0.69;   % 0.71                0.69256 , participants look at the E option first

%% Updated Fixations from 25.08.2025
if OV == 1
    FixR = rand() > (0.70*(Vr>Vl) + 0.67*(Vr<Vl) + 0.67 *(Vr==Vl) ); 
    FixDur = [round(lognrnd(5.26, 0.59, 1, 1)); round(lognrnd(5.37, 0.42, pre_Nfix-1, 1))];

elseif OV == 2
    FixR = rand() > ( 0.68*(Vr>Vl) + 0.68 *(Vr<Vl) + 0.65 *(Vr==Vl) );
    FixDur = [round(lognrnd(5.30, 0.63, 1, 1)); round(lognrnd(5.30, 0.44, pre_Nfix-1, 1))];

elseif OV == 3
    FixR = rand() > ( 0.70*(Vr>Vl) + 0.70*(Vr<Vl) + 0.71*(Vr==Vl) );
    FixDur = [round(lognrnd(5.31, 0.58, 1, 1)); round(lognrnd(5.34, 0.46, pre_Nfix-1, 1))];
end

% this can be chosen instead of the n-bcak fixation mechanism 
% if FixR
%     FixR(2:pre_Nfix+1) = repmat([1 0]',pre_Nfix/2,1);
% else
%     FixR(2:pre_Nfix+1) = repmat([0 1]',pre_Nfix/2,1);
% end

% empirical P(fixate E) as a function of fix index from end,
% that is, from point of a decision
% for ES phase data

pE_from_end = [ ...
    0.47 ... % at choice
    0.53 ... % 1 before choice
    0.61 ... % 2 fix beforee choice
    0.62 ... % 3 fixations before choice  (on average, there's ca 3.2 fix)
    0.59];   % 4 fixations before

max_back = numel(pE_from_end);

FixR = nan(pre_Nfix,1);  % 0 = left, 1 = right

for k = 1:pre_Nfix
    % index counted backwards from decision
    idx_back = pre_Nfix - k;   % 0 = last fixation, 1 = one before, etc.

    if idx_back <= max_back-1
        p_fixE = pE_from_end(idx_back + 1);
    else
        p_fixE = pE_from_end(end); 
    end

    if rand < p_fixE
        FixR(k) = 0;   % fixate left
    else
        FixR(k) = 1;   % fixate right
    end
end

%-----------------------------------------------------------------------

k_fix    = 0;

%% pre-allocation
dvPeriod = cell(pre_Nfix,1); % to save multiple dv for each fixation
dv       = zeros(pre_Nfix,1)+spbias;  % + spbias, can be added (select depending on model)

%% evidence accumulation for each fixation
while k_fix< pre_Nfix
    
    k_fix = k_fix+1;

    % for ESaDDM
%     if FixR(k_fix) == 1
%         dv(k_fix) = beta2.*(Vr-Vl) - beta3.*(Vl-Vr); 
%     else
%         dv(k_fix) = beta4.*(Vr-Vl) - beta2.*(Vl-Vr); 
%     end

%aDDM
    if FixR(k_fix) == 1
        dv(k_fix) = beta2.*(Vr-Vl)-beta3.*(Vl-Vr); 
    else
        dv(k_fix) = beta3.*(Vr-Vl)-beta2.*(Vl-Vr); 
    end



    err = normrnd(0,noisy,1,FixDur(k_fix)); 
    
    dvPeriod{k_fix,1} = dv(k_fix)+err; 
    
    Fix.A{k_fix,1} = repmat(FixR(k_fix),1,FixDur(k_fix)); 
    
end

%% sum up dv

rv   = [dvPeriod{1:pre_Nfix}];
sumdvALL = cumsum(rv)+spbias;      % spbias +  with the z bias  

% identify the dv crossing threshold
ID = find(abs(sumdvALL)>=a);

if isempty(ID)
    display('Not enough evidence')
    E                 =nan; % evidence at the time point crossing threshold
    tempEyeData.DwellR        = nan;
    tempEyeData.DwellL        = nan;
    tempEyeData.Dwelltotal    = nan;
    tempEyeData.Nfix   = nan;
    tempEyeData.FixLocFirst    = nan;
    tempEyeData.FixLocLast     = nan;
    tempEyeData.FixLocFirstCorr    = nan;
    tempEyeData.FixLocLastCorr     = nan;
    tempEyeData.FirstFixDur      = nan;
    tempEyeData.MiddleFixDur     = nan;
    tempEyeData.FinalFixDur      = nan;
    tempEyeData.eachMiddleFixDur = nan;
    tempEyeData.DwellDiff = nan;
    FixAaLL = nan; sumdvALL = nan; FixDur = nan;
    Choice = nan;
    RT     = nan;
else
    %% organize simulated data into the same array
    FixAaLL      = [Fix.A{1:k_fix}]; 
    FixationInfo =  FixAaLL(1:ID(1));
    
    %% summarize data
    E                 = sumdvALL(1:ID(1)); % evidence at the time point crossing threshold
    tempEyeData.DwellR        = sum(FixationInfo(1,:) == 1);
    tempEyeData.DwellL        = sum(FixationInfo(1,:) == 0);
    tempEyeData.DwellDiff     = tempEyeData.DwellR-tempEyeData.DwellL;
    tempEyeData.Dwelltotal    = tempEyeData.DwellR +tempEyeData.DwellL;
    tempEyeData.Nfix          = sum(diff(FixationInfo(1,:))~=0)+1; % number of real fixation (only count switching between A and B)
    
    if tempEyeData.Nfix > 1
        tempEyeData.FixLocFirst   = FixR(1);
        tempEyeData.FixLocLast    = FixationInfo(1,end-1);
        tempEyeData.FixLocFirstCorr   = FixR(1) == (Vr>Vl);
        tempEyeData.FixLocLastCorr    = tempEyeData.FixLocLast == (Vr>Vl);
    else
        tempEyeData.FixLocFirst   = FixR(1);
        tempEyeData.FixLocLast    = NaN;
        tempEyeData.FixLocFirstCorr   = FixR(1) == (Vr>Vl);
        tempEyeData.FixLocLastCorr    = NaN;
    end
    
    
    if tempEyeData.Nfix > 2
        SwitchPoints              = find(diff(FixationInfo(1,:)));
        
        tempEyeData.FirstFixDur      = SwitchPoints(1);
        tempEyeData.FinalFixDur      = tempEyeData.Dwelltotal - SwitchPoints(end); % final prefixation duration + (Dwell - total prefixation duration)
        tempEyeData.MiddleFixDur     = tempEyeData.Dwelltotal - tempEyeData.FirstFixDur - tempEyeData.FinalFixDur; % Nfix - final fixation - first fixation duration
        tempEyeData.eachMiddleFixDur = tempEyeData.MiddleFixDur./(tempEyeData.Nfix-2);
        
        
    elseif tempEyeData.Nfix == 2
        SwitchPoints              = find(diff(FixationInfo(1,:)));
        
        tempEyeData.FirstFixDur      = SwitchPoints(1);
        tempEyeData.MiddleFixDur     = nan;
        tempEyeData.FinalFixDur      = tempEyeData.Dwelltotal - SwitchPoints(end);
        tempEyeData.eachMiddleFixDur = nan;
    elseif tempEyeData.Nfix == 1
        tempEyeData.FirstFixDur      = tempEyeData.Dwelltotal;
        tempEyeData.MiddleFixDur     = nan;
        tempEyeData.FinalFixDur      = nan;
        tempEyeData.eachMiddleFixDur = nan;
    end
    
    
    Choice = E(end)>0; % 1: Right; 0: Left
    RT     = length(E)+(ndt*1000); %*1000
    %RT     = length(E)+ ndt; %*1000
end