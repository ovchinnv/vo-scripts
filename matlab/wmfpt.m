% compute MFPT from F and D
%
fefile='wfe.mat';
diffile='wdiff.mat';

load(fefile);
load(diffile);
fe1=fe;
% to use fe computed from cfe2i (should make no significant difference) :
load 'fe-rm5b-rm5a.mat' ;

%plot(mean(fe1,1)); hold on;
%plot(mean(fe,1));
%return

% test using average D :
%Diff(:)=mean(Diff);
%Diffe(:)=mean(Diffe);

kBT=kboltz*Temp;

[nsamp,np]=size(fe);
ferr=var(fe,1)/nsamp; ferr(1:np)=0;
fe=mean(fe,1);
%%%%%%%%% compute MFPT integrals
%
% first compute inner integral
fe=fe-fe(1) ;
%pdfc=0.5 * ( exp ( -fe(1:end-1)/kBT ) + exp ( -fe(2:end)/kBT ) ) ;
pdfc = exp ( -fe(1:end-1)/kBT ) ;
pdfce = ( pdfc/kBT ).^2 .* ferr(1:end-1) ; % variance
%pdfc = exp ( -fe(2:end)/kBT ) ;

dds=diff(alpha) ;
%
i1(1)=0; % first point of inner integral
i1e(1)=0;
for i=1:np-1
 i1(i+1)  = i1(i)  + pdfc(i)  .* dds(i)  ;
 i1e(i+1) = i1e(i) + pdfce(i) .* dds(i)  ; % error
end

%i1c=0.5*(i1(1:end-1) + i1(2:end)) ; % shift to center
i1c=i1(1:end-1) ;
i1ce=i1e(1:end-1) ; % error
%i1c=i1(2:end) ;
% second integrand :
%pdfcd=0.5 * ( exp ( fe(1:end-1)/kBT )./Diff(1:end-1) + exp ( fe(2:end)/kBT )./Diff(2:end) ) ;
pdfcd = exp ( fe(1:end-1)/kBT )./ Diff(1:end-1) ;
pdfcde = ( pdfcd/kBT ).^2 .* ferr(1:end-1) + (pdfcd./Diff(1:end-1)).^2 .* Diffe(1:end-1) ; %error
%pdfcd = exp ( fe(2:end)/kBT )./ Diff(2:end) ;

i2 = i1c.*pdfcd.*dds ; % summing this up from index i to j gives the required mfpt ; 
i2e = ( i1ce .* pdfcd.^2 + i1c.^2 .* pdfcde ) .*dds ;

% plot MFPT from a point to the end :
ifrom=1 ;
ito=31 ;
mfpt=[cumsum(i2(ifrom:end))];
mfpte=[cumsum(i2e(ifrom:end))];
%
figure(1) ; 
%clf ;
hold on; box on;
%plot([ifrom+1:ito],mfpt(1:ito-ifrom),'r.-');
alpha2=( alpha(ifrom+1:ito) - alpha(ifrom) ) / ( alpha(ito) - alpha(ifrom) )
%errorbar(alpha(ifrom+1:ito),mfpt(1:ito-ifrom),sqrt(mfpte(1:ito-ifrom)),'k.--','markersize',12);
%plot(alpha2,mfpt(1:ito-ifrom),'r-','markersize',12);
errorbar(alpha2,mfpt(1:ito-ifrom),sqrt(mfpte(1:ito-ifrom)),'k.--','markersize',12);

set(gca,'yscale','log')

grid on ;
%set(gca,'ytick',[0 : 1e-6 : 1e-5]);
ylabel('\it MFPT(s)', 'fontsize',14);
%xlabel('\it \alpha', 'fontsize',14);
xlabel('\it Reaction coordinate (\alpha)', 'fontsize',14);
%
% label states
text(0.0,2,'\it \bf M5B');
text(0.9,2,'\it \bf M5A');
text(0.5,5.5,'\it \bf FE barrier');
title('\it \bf M5A bound to rab11a in two conformations (M5B - M5A)');
%
xlim([0 1]);
%ylim([0 1e-5]);
set(gcf, 'paperpositionmode', 'auto');


%print(gcf, '-depsc2','wmfpt') ;

