function [eegbb, trig01, trigtype, eegmb] = EphysEpochFxn(path,trigch,dur)

%% extraction of lfp, csd, mua
% this function takes a wavelet file path and epochs the MUA and LFP 
% edits by Chase Mackey 2022


% load data
%path = ('C:\Users\cmackey\Documents\MATLAB\1-gt006000010@os_eye06_30');
load(path);
fname = path;
cntm = wraw.cntm; %MUA
%cntb = wraw.cntb; %LFP
eegbb=1;

%% INPUTS (filtering, epoch timing etc.)%
timeframe_baseline              = [-dur 0];
epoch_tframe            = [-dur dur];  
newadrate =500;

%

%% organize data based on input trigger "trigch"

% set timeframe
time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);


% aud - 1, vis (open eyes) - 3, sacc on 4, sacc off 5
trig0=trig.anatrig{trigch};


%triggers if there are multiple, trigch is the one we're selecting


if trigch==1  % aud stim is stored as cell type SOMETIMES so we need a lot of ifs
    
    if iscell(trig.ttype)
        present_trigidx=trig.ttype{1,trigch}(:,:)>0; % places where there was a trigger
        trigtype=trig.ttype{1,trigch}(present_trigidx); % aud stim are in cells usually
    else
        present_trigidx=trig.ttype(:,:)>0;
        trigtype=trig.ttype(present_trigidx); 
    end
    
elseif trigch==2
    present_trigidx=trig.ttype(:,:)>0;
    trigtype=trig.ttype(present_trigidx); % vis stim are in double ... should be fixed
    
elseif trigch==3 % led stim when eyes are open, gotta do something special here
    
    trig0_all=trig.anatrig{2}; %ALL (eyes open & close) trigs
    for i=1:length(trig0)
       trig_mutual(i)= find(trig0_all==trig0(i)); %gives you the ALL trig number that is in trig0
    end
    trigsorttmp=trig.triglength{1,2};
    trigtype(:,1)=trigsorttmp(trig_mutual); % std/dev when eyes are open % vis stim are in double ... should be fixed
 
    
elseif trigch==4 % sacc on
    trigtype=ones(size(trig.anatrig{1,trigch}(:)));% just sticking ones here because we don't have diff sacc types
    
elseif trigch==5 % sacc off
    trigtype=ones(size(trig.anatrig{1,trigch}(:))); % just sticking ones here because we don't have diff sacc types
end


%% sorting by trigger type
% uniquetrigs = unique(trigtype);

% if length(uniquetrigs)>1 && length(uniquetrigs)<3
%     
%     stds = trigtype==uniquetrigs(1)';
%      %devs = trigtype==uniquetrigs(2);
%     trig0 = trig0(stds(1:length(trig0))==1); % trigtype/trig0 mismatch?
%     
% end

% epoch to the std stimulus, or the only trig if there's only one
for trigredxct=1:length(trig0)
    trig01(trigredxct)    = round(trig0(trigredxct)./(trig.adrate/newadrate));
end

%% FILTERING 
%filter the cnt in 1-100Hz to make nice CSD picts
%to bandpass filter in the 1.2-2.5Hz range the MUA (cntm)
% n = 2;
% Wn = [1 100]/(newadrate/2); %the first 2 numbers are the freq range you want to pass, which is =/- 20% of the reprate
% [b,a] = butter(n,Wn);
% 
% % preallocate for filtered data
% % cntc_ff = zeros(size(cntm));

% % cnte_ff = zeros(size(cntm));
% cntb_ff = zeros(size(cntm));
% 
% i=1;
%  for i=1:size(cntm,1) %this is for the electrode chs
% % %  cntc_ff(i,:)=filtfilt(b,a,cntc(i,:)); %creates a new variable which contains filtered data
% % %  cnte_ff(i,:)=filtfilt(b,a,cnte(i,:));
% %  %cntb_ff(i,:)=filtfilt(b,a,cntb(i,:));
%  end

% epoch times
x1 = round(epoch_tframe(1)*(newadrate/1000));
x2 = round(epoch_tframe(2)*(newadrate/1000));

%% sometimes have too many triggers 

ExpectedLastTrialEnd = trig01(length(trig01))+x2; % last trial based on trigs
SkipThese = zeros(1,length(trig01)); % Start by skipping none

while ExpectedLastTrialEnd > length(cntm) % too many triggers, 
    if exist('lastreal')
        trig01(1,lastreal) = NaN; % so mark last trial with NaN
        SkipThese = isnan(trig01(1,:)); % get idx of last trial
    else
        trig01(1,length(trig01)) = NaN; % so mark last trial with NaN
        SkipThese = isnan(trig01(1,:)); % get idx of last trial
    end
    
    lastreal = length(trig01(1,SkipThese==0)); % last good trig time
    ExpectedLastTrialEnd = trig01(lastreal)+x2; % new expt end time

    if ExpectedLastTrialEnd <= length(cntm) % break out 
        break
    end


end


%% epoching
i=1;
for i=1:length(trig01(1,SkipThese==0)) %skipping extra trigs when epoching
    
 eegm(:,i,:)=cntm(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
%  eegc(:,i,:)=cntc_ff(:,trig01(i)+x1:trig01(i)+x2);%epoching the csd
%  eege(:,i,:)=cnte_ff(:,trig01(i)+x1:trig01(i)+x2);% same for lfp
 %eegb(:,i,:)=cntb_ff(:,trig01(i)+x1:trig01(i)+x2);% same for bipolar lfp
 
end


%% baseline correct
chct=1; % channels
trct=1; % trials

for chct=1:size(eegm,1)
 for trct=1:size(eegm,2)
     eegmb(chct,trct,:)   = squeeze(eegm(chct,trct,:))-squeeze(mean(eegm(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     %eegcb(chct,trct,:)   = squeeze(eegc(chct,trct,:))-squeeze(mean(eegc(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     %eegeb(chct,trct,:)   = squeeze(eege(chct,trct,:))-squeeze(mean(eege(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
    % eegbb(chct,trct,:)   = squeeze(eegb(chct,trct,:))-squeeze(mean(eegb(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
 end
end




end