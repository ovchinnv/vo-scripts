% read text file with standard deviation of the force and compute a window
% index that is expected to reduce overall force error
%

s=load('dfe.dat');
amp=2; % amplification exponent : 1 would reproduce s
sn=cumsum(s.^amp); sn=sn/sn(end);

% treat sn as a cdf and use it to randomly sample from s, interpreted as a pdf
nsamp=1000;
nsamp=10;
r=rand(1,nsamp);
iwins=[];
for i=1:nsamp
 iwin=find(sn>r(i),1);
 iwins=[iwins iwin];
end

%write windows to file :
save -ascii 'nextwin.dat' iwins

return;

close ;
[h,x]=hist(iwins,[1:length(s)]);

plot(s/sum(s),'k-x'); hold on
plot(x,h/sum(h)/(x(2)-x(1)),'r:');

