% time series bootstrap following Kunsch 1989
% the error that I get is higher than that from the bootstrap matlab package
%
% to do :  rewrite as a function

x=pbe1(:);

m=100 ; % provided sample decorrelation length

N=length(x);
n=N-m+1;

C=1;
l=C * floor(n^0.333) ; % sampling block length

l=m

% number of blocks of length l
k=floor(n/l); 
% recompute n such that k * l = n
n=k*l;

nsample=999; % number of resampling trials

% set seed for reproducibility
rand('seed', pi) 

%sample block offsets from uniform distribution
S = floor ( rand(k , nsample) * (n-l+1) ) ;

% sample the time series according to S : 
f=zeros(1,n);
fmean=zeros(1,nsample);
%
for i = 1 : nsample
 q=1;
 for j = 1 : k
  ind = S ( k , i);
  for t = 1 : l
% evaluate function
% 1 - get series
   xx=x ( ind + t : ind + t + m - 1 );
% 2 - evaluation
%     create new, presumably uncorrelated, time series with observable of interest
   f(q)=mean(xx) ; % mean of the series (could be other statistics, e.g. standard deviation)
   q=q+1;
  end
 end
%
 fmean(i)=mean ( f ) ; % mean of resampled
%
end

mean(fmean)
std(fmean)

