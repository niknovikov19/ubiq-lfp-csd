function [eegcb, eegbb, trig01, trigtype, eegmb,epoch_tframe,std_freq,eege] = EphysExtractFxn(path,trigch)

%% extraction of lfp, csd, mua
% this function takes a file path and gets LFP and MUA data out of the file
% edits by Chase Mackey 2022


% load data
%path = ('C:\Users\cmackey\Documents\MATLAB\1-gt006000010@os_eye06_30');
load(path);
fname = path;


%% INPUTS (filtering, epoch timing etc.)%
timeframe_baseline              = [-30 0];
epoch_tframe            = [-300 300];  
newadrate =1000;
filtere = [0.5 300];%LFP
filteru = [300 5000];%MUA
filtertype=1;
fsize = 6;
xlabelres=5;
%

%% organize data based on input trigger "trigch"

% set timeframe
time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);

% extract continuous ephys data and triggers
[~, cnte, cntm, cntc, ~, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);

%[trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig, params);

% noise - 1, led (open eyes) - 3, sacc on 4, sacc off 5
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


for trigredxct=1:length(trig0)
 trig01(trigredxct)    = round(trig0(trigredxct)./(trig.adrate/newadrate));
end

% if trig type doesn't match trig0 for some reason (happens..need to figure
% it out
% Find the indices in trig01 that don't exceed the size of cntm
valid_indices = trig01 <= size(cntm, 2);

% Trim trig01 to keep only the valid indices
trig01 = trig01(valid_indices);

% tone frequency HARD CODED FOR NEUROSCAN FILES
if length(params.filedat)>15
    std_freq = params.filedat(16);
else
    disp('tone freq not found, is this not a neuroscan file?')
    std_freq = 0;
end
%% 
%filter the cnt in 1-100Hz to make nice CSD picts
%to bandpass filter in the 1.2-2.5Hz range the MUA (cntm)
n = 2;
Wn = [1 100]/(newadrate/2); %the first 2 numbers are the freq range you want to pass, which is =/- 20% of the reprate
[b,a] = butter(n,Wn);

% preallocate for filtered data
cntc_ff = zeros(size(cnte));
cntm_ff = zeros(size(cnte));
cnte_ff = zeros(size(cnte));
cntb_ff = zeros(size(cnte));

i=1;
% for i=1:size(cntc,1) %this is for the electrode chs
%      cntc_ff(i,:)=filtfilt(b,a,cntc(i,:)); %creates a new variable which contains filtered data
%      cntm_ff(i,:)=filtfilt(b,a,cntm(i,:));
%      cnte_ff(i,:)=filtfilt(b,a,cnte(i,:));
%      cntb_ff(i,:)=filtfilt(b,a,cntb(i,:));
% end

%epoching
x1 = round(epoch_tframe(1)*(newadrate/1000));
x2 = round(epoch_tframe(2)*(newadrate/1000));

i=1;
for i=1:length(trig01)
    
 eegm(:,i,:)=cntm(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
 eegc(:,i,:)=cntc(:,trig01(i)+x1:trig01(i)+x2);%epoching the csd
 eege(:,i,:)=cnte(:,trig01(i)+x1:trig01(i)+x2);% same for lfp
 eegb(:,i,:)=cntb(:,trig01(i)+x1:trig01(i)+x2);% same for bipolar lfp
 
end

% for i=1:length(trig01)
%     
%  eegm(:,i,:)=cntm(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
%  eegc(:,i,:)=cntc(:,trig01(i)+x1:trig01(i)+x2);%epoching the csd
%  eege(:,i,:)=cnte(:,trig01(i)+x1:trig01(i)+x2);% same for lfp
%  eegb(:,i,:)=cntb(:,trig01(i)+x1:trig01(i)+x2);% same for bipolar lfp
%  
% end


%% baseline correct
chct=1; % channels
trct=1; % trials

for chct=1:size(eegm,1)
 for trct=1:size(eegm,2)
     eegmb(chct,trct,:)   = squeeze(eegm(chct,trct,:))-squeeze(mean(eegm(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegcb(chct,trct,:)   = squeeze(eegc(chct,trct,:))-squeeze(mean(eegc(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegeb(chct,trct,:)   = squeeze(eege(chct,trct,:))-squeeze(mean(eege(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegbb(chct,trct,:)   = squeeze(eegb(chct,trct,:))-squeeze(mean(eegb(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
 end
end




end