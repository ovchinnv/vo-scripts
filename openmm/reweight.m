% reweight time series to lowest temperature
% compute the log of partition functions first, using a restart file
% read time series and compute weighting factor
%
T0=305 ;% target temperature

% replicas :
brep=1;
erep=15;

% over all replicas : 
for irep=brep:erep
fprintf([' Processing replica ', num2str(irep),'\n']);
%
% iterations :
brun=1;
erun=100;
%
basename=['../pmfsep1_',num2str(irep),'.temp_'];
%
restart_file = [basename,num2str(erun),'.txt'];
% optional file to remove earlier statistics ; comment to disable
%restart_subtract = 'hhh_4.0_05.temp15.txt';
gdata;

close all;

% interpolate average using the method of Zhang&Ma10:
n=length(bet) ;
db=50; % number of interpolation bins on either side of central bin
for i=1:n
 ib=max(1, i-db);
 ie=min(i+db,n);
 x=bet (ib:ie)';
 y=eavg(ib:ie)';
 w=wgt (ib:ie)';
% determine fit
 m=length(x);
 ws=sum(w);
 wx=w.*x;
 xs=sum(wx) ; 
 ys=w*y' ;
 xy=wx*y' ;
 x2s=wx*x' ;
 d = ws*x2s-xs^2 ;
 a = (x2s*ys-xs*xy)/d;
 b = (ws*xy-xs*ys)/d;
% evaluate :
 eint(i)=a+b*bet(i) ;
end
% integrate the interpolated average energy to obtain log ratio of partition functions Z(beta)/Z(beta(end))
% note that the highest beta corresponds to the lowest temperature
% use trapezoid rule
logZ(n)=0; % initialize
for i=n-1:-1:1
 logZ(i) = logZ(i+1) + 0.5 * (-1) * ( eint(i) + eint(i+1) ) * ( bet(i) - bet(i+1) );
end
if (0)
% plot logZ :
figure ; hold on; box on;
plot(bet, logZ,'k', 'linewidth',2);
ylabel('\it log Z(\beta)/Z(\beta_0)', 'fontsize',14);
xlabel('\it \beta(kcal/mol)^{-1}', 'fontsize',14);
end
%
% now read time series files and compute weights ; add an extra column, and write out a new file
%

files={};
for i=brun:erun
 files=[files {[basename,num2str(i),'.series.txt']}];
end

for i=1:length(files)
 file=char(files(i));
 if (i==1)
  d=load(file);
 else
  d=[d; load(file)];
 end 
end

% instantaneous values :
istep=d(:,1);
itemp=d(:,2);
iener=d(:,3);
% compute inst. beta :
ibet=1./(itemp*kb);
% compute logZ/logZ(0) by interpolation :
ilogZ=interp1(bet,logZ,ibet,'linear','extrap');
% target beta : 
bet0=1./T0/kb;
logZ0=interp1(bet,logZ,bet0,'linear','extrap');

ilogw=(ibet-bet0).*iener + ilogZ - logZ0 ;

ilogw=min(ilogw,0);

%plot(ilogw);
%plot(exp(ilogw))
% write series file with log weights
fname=[basename,'series.txt'];
%d=[istep, ibet, itemp, iener, ibet-bet(n), ilogZ, ilogw, exp(ilogw)];
d=[istep, itemp, iener, ilogw];
save(fname,'-ascii', 'd');

end % irep


