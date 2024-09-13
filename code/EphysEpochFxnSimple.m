function [lfp] = EphysEpochFxnSimple(cnt,trig01,epoch_tframe,newadrate)
%% simple epoch fxn
% inputs: continuous data, triggers in samples, desired timeframe in ms

x1 = round(epoch_tframe(1)*(newadrate/1000));
x2 = round(epoch_tframe(2)*(newadrate/1000));
i=1;
for i=1:length(trig01)
    
    lfp(:,i,:)=cnt(:,trig01(i)+x1:trig01(i)+x2);% same for lfp

end