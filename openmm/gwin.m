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
wins=[];
for i=1:nsamp
 win=find(sn>r(i),1)+iwin-1; % offset to zero
 wins=[wins win];
end

%write windows to file :
save -ascii 'nextwin.dat' wins

return;

close ;
[h,x]=hist(wins,[1:length(s)]);

plot(s/sum(s),'k-x'); hold on
plot(x,h/sum(h)/(x(2)-x(1)),'r:');

