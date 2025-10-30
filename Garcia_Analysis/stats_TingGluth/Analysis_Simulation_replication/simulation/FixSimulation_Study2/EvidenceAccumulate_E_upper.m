% - Main script used to perform evidence accumulation for each trial
% 2022/12: add output: DwellFirst, DwellMiddle, DwellFinal, FixLocFirstCorr
% - modified from Dr Chih-Chung Ting's script
% IMPORTANT: This code assumes response coding with the upper bound being
% the E option (for analysisi purposes E is always the left option; so when
% chose_left == 1, E was chosen


function [Choice, RT, E, tempEyeData, FixAaLL, sumdvALL, FixDur]= EvidenceAccumulate_E_upper(OV,Vl,Vr,params)

% parameters
a     = params(1).*0.5;        % boundary seperation 
ndt   = params(2);             % non-decision time
%z     = params(3);            % starting point
beta2 = params(3)./1000;       % intercept
beta3 = params(4)./1000;       % attention or for ddm Value_difference
%beta3 = params(5)./1000;       % inattention
%beta4 = params(7)./1000;
%beta5 = params(8)./1000;

noisy = 0.03;                  %params(6);       %
%theta = beta3./beta2;
% z 
%spbias = (2*z - 1) * a;       % if no bias, then z=0.5, and spbias simply 0, but I am unsure if this is the correct way of adding z


%% pre-determine fixations
pre_Nfix = 70; % assume 70 fixations before decision
FixR     = rand()> 0.69;   % 0.71                0.69256 , participants look at the E option first

% Updated Fixations from 25.08.2025
% because participants look at the E (left) option ca 70% of the time, the
% probabiltiy of fixating right is always pretty low
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


if FixR
    FixR(2:pre_Nfix+1) = repmat([0 1]',pre_Nfix/2,1);
else
    FixR(2:pre_Nfix+1) = repmat([1 0]',pre_Nfix/2,1);
end


%% initializing the fixation process
k_fix    = 0;

% pre-allocations
dvPeriod = cell(pre_Nfix,1); % to save multiple dv for each fixation
dv       = zeros(pre_Nfix,1);  

% evidence accumulation for each fixation
while k_fix< pre_Nfix
    
    k_fix = k_fix+1; % Fixation from fixation

        if FixR(k_fix) == 1
            dv(k_fix) = beta2.*Vr-beta3.*Vl; % for each ms. In HDDM, the smallest unit is s.
        else
            dv(k_fix) = beta3.*Vr-beta2.*Vl; % for each ms. In HDDM, the smallest unit is s.
        end


% % ESaDDM function:
%         if FixR(k_fix) == 1
%             dv(k_fix) = beta2.*(Vr-Vl) - beta3.*(Vl-Vr); 
%         else
%             dv(k_fix) = beta4.*(Vr-Vl) - beta2.*(Vl-Vr); 
%         end

    
    err = normrnd(0,noisy,1,FixDur(k_fix)); % Sampled errors for N (FixDur(k_fix)) time points.
    
    dvPeriod{k_fix,1} = dv(k_fix)+err; % dv for each time point within a fixation (the combination of dv and error).
    
    Fix.A{k_fix,1} = repmat(FixR(k_fix),1,FixDur(k_fix)); % label each time point as fixating at the better option A.
    
end

%% sum up dv

rv   = [dvPeriod{1:pre_Nfix}];
sumdvALL = cumsum(rv);      % spbias +  with the z bias - CCT version



% identify the dv crossing threshold
ID = find(abs(sumdvALL)>=a);

if isempty(ID)
    display('Not enought evidence')
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
    FixAaLL      = [Fix.A{1:k_fix}]; % put all lables (fixating at the better or worse option) into the same array
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
end


