% also, output average energies vs T for all replicas in one file
%
if ~exist('styles')
 styles={'r-s','g-x','b-o','m-*','c-v','k-^','r--','g--','b--','m--','c--','k--'};
end

%
kboltz=1.987191e-3 ;
% replicas :
brep=1;
erep=15;
% restart file index
erun=100;

% over all replicas :
for irep=brep:erep
fprintf([' Processing replica ', num2str(irep),'\n']);
%
% read restart file :
basename=['../pmfsep1_',num2str(irep),'.temp_'];
restart_file = [basename,num2str(erun),'.txt'];
% optional file to remove earlier statistics ; comment to disable
%restart_subtract = 'hhh_4.0_05.temp15.txt';
gdata;

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

eave(:,irep+1)=eint;

end % irep
%
eave(:,1)=bet ; % temperature

% write energies

outname=['pmfsep1_ene_temp',num2str(erun),'.dat'];
save(outname,'-ascii', 'eave');

% plot landscape
temp=1./bet/kboltz;
%pcolor(brep:erep,temp,eave(:,brep+1:erep+1)) ; shading interp ;

% plot several temperature curves
temps=[300 : 25 : 375];
leg={};
%
close all;
fig=figure; hold on; box on;
for k=1:length(temps)
 t=temps(k)
 for irep=brep:erep
  ethis(irep)=interp1(temp,eave(:,irep+1),t, 'linear','extrap');
 end
%
 plot(brep:erep,ethis-ethis(1), char(styles(k)))
%
leg=[ leg {['T=',num2str(t)]}];
end

% plot average energy across all temperatures (note that this is actually an average over all beta; this is important because the temperature
% and beta pdfs are _different_)
eave_bet=mean(eave(:,brep+1:erep+1),1);
plot(brep:erep,eave_bet-eave_bet(1),'k*-','linewidth',2)
leg=[leg {['Average over \beta']}];

legend(leg,4)

ylabel('\it <E>(kcal/mol)', 'fontsize',14);
xlabel('\it Window #','fontsize',14);
