X = load('RelPowAreaA1.mat');
X = X.RelPow;
X(5) = [];

%X = load('RelPowAreaA1_CSD.mat');
%X = X.results;

Y = load('RelPowAreaA1.mat');
l4_chans = [Y.RelPow.L4chan];
l4_chans(5) = [];

nchan = size(X(1).relpow, 1);
nfreq = size(X(1).relpow, 2);
nsess = length(X);

nchan_max = 35;
ch0 = 17;

ch1 = ch0 - round(nchan/2);
ch2 = ch1 + nchan - 1;

Q = NaN * ones(nchan_max, nfreq, nsess);

for n = 1 : nsess   
    d = ch0 - l4_chans(n);
    Q(ch1+d : ch2+d, :, n) = X(n).relpow;    
end

W = nanmean(Q, 3);
figure(); hold on;
imagesc(W);
plot([1, nfreq], [ch0, ch0], 'k--');
set(gca, 'YDir', 'reverse');
colorbar;