function [Metrics] = ModelMetrics(x1,y1,lim,w)
% This function to calculate Jeffrey Divergence, Kullback-Liebler
% Divergence and calculate Match Distance between P and Q.
% All as difined at [1]
%
%   [1] Rubner, Yossi, Carlo Tomasi, and Leonidas J. Guibas. "The earth mover's 
%   distance as a metric for image retrieval." International journal of 
%   computer vision 40.2 (2000): 99-121.
%
%
%   x - Time domain signal from patient;
%   y - Time domain signal from model

%Getting P(x) and Q(x)



Metrics=struct;
hist1=histogram(x1,'BinLimits',lim,'BinWidth',w,'Normalization','probability');
h=hist1.BinCounts;
hold on
hist2=histogram(y1,'BinLimits',lim,'BinWidth',w,'Normalization','probability');
k=hist2.BinCounts;

if (length(h)> length(k))
    h=h(1:length(k));
end

if (length(k)> length(h))
    k=k(1:length(h));
end

h=h./sum(h); %normalization
k=k./sum(k);

%% Intersection Distance




num=[];
for i=1:length(h)
    num(i)= min(h(i),k(i));
end
Metrics.dI=(sum(num)/sum(h));

%% L-1 distance (Manhattan)
H=cumsum(h);
K=cumsum(k);
Metrics.dL1=sum(abs(H-K));



%% Jensen-Shannon Divergence (link: encr.pw/8gbtV)
M=0.5*(h+k);

dKL1=(h.*log2(h./M)); %% (P||M)
dKL1(isnan(dKL1))=0; % 0/0 problem
DKL1=sum(dKL1);

dKL2=(k.*log2(k./M)); %% (Q||M)
dKL2(isnan(dKL2))=0; % 0/0 problem
DKL2=sum(dKL2);

JSD=0.5*DKL1+0.5*DKL2;

Metrics.JSD=JSD;



%% Kullback-Liebler Divergence
k(k==0) = 1e6; % division by zero problem

dKL=(h.*log2(h./k));
dKL(isnan(dKL))=0; %0/0 problem
Metrics.dKL=sum(dKL);



%% Stat Moments
Metrics.var=[var(x1) var(y1)];

Metrics.Kurtosis=[kurtosis(x1) kurtosis(y1)];
Metrics.Skewness=[skewness(x1) skewness(y1)];

%centroid
[yip,xip] = histcounts(x1);
[yiq,xiq] = histcounts(y1);
yip(end+1)=yip(end);
yiq(end+1)=yiq(end);
Metrics.centroid_P=[sum(xip.*yip)/sum(yip) sum(yip.*xip)/sum(xip)];
Metrics.centroid_Q=[sum(xiq.*yiq)/sum(yiq) sum(yiq.*xiq)/sum(xiq)];



end

