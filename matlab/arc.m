%
close all;
addpath /home/taly/scripts/matlab;

if (~exist('read'))
 read=1;
end

if (read==1)
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'arcl2.dat' 'arcl3.dat' 'arcl4.dat'  ...
 }
%
clear s;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  s=load(fname);
 else
  s=[s; load(fname)];
 end 
end
read=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[niter,nstring]=size(s);

ind=[1:niter];
d=1;
ssum = sum(s(:,2:end),2);
ssum = smooth2(ind, ssum, d);
%
%
plot(ind, ssum );

box on;
ylabel('Replica index');
xlabel('iteration');

set(gcf, 'paperpositionmode', 'auto');

%print(gcf, '-depsc2', 'slen.eps');
%print(gcf, '-djpeg100', 'slen.jpg');



