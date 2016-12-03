% plot all rmsf in one figure 
% 
close all;

names={ '3bnc60t' '3bnc60at', '3bnc60glt',};
flags={'-ca-hc', '-ca-lc' };
flags={'-cg-hc', '-cg-lc' };
clrs={'r','g','b'};
aligns={'malign.mat' 'malign-lc.mat'};
xlims=[125,110];

lw=1;
%
for jj=1:2
  leg={};
  hleg=[];
  flag=char(flags(jj));
  align=char(aligns(jj));
  load(align);
  armsf=zeros(size(maind));
% indicate domains
  figure('position',[100 300 900 300]); hold on; box on;
  showdom;
% plot RMSFs
  for ii=1:3
    armsf(:)=NaN; % do not forget to reinitialize
    name=char(names(ii));
    clr=char(clrs(ii));
    load([name,'-rmsf',flag,'.mat']);
    n=length(flucall);
    resnum=[1:n];

    addpath('~/scripts/matlab');
    d=3;
    fs =smooth2(resnum,flucall,d);
    fsp=smooth2(resnum,flucall+flucstd,d);
    fsm=smooth2(resnum,flucall-flucstd,d);
%
% populate armsf
    for k=1:length(armsf)
     i3=maind(k,ii);
     if (i3>0 & i3<=n)
      armsf(k,1)=fs(i3);
      armsf(k,2)=fsp(i3);
      armsf(k,3)=fsm(i3);
     end
    end

%    h=plot(resnum,fs,clr,'linewidth',lw) ;hold on ; box on ;
%    plot(resnum,fsp,[clr,'--'],'linewidth',1) ;hold on ; box on ;
%    plot(resnum,fsm,[clr,'--'],'linewidth',1) ;hold on ; box on ;
  h=plot(armsf(:,1), clr,'linewidth',lw) ;hold on ; box on ;
    plot(armsf(:,2),[clr,'--'],'linewidth',1) ;hold on ; box on ;
    plot(armsf(:,3),[clr,'--'],'linewidth',1) ;hold on ; box on ;
    leg = [ leg {[name,flag]} ];
    hleg=[hleg h]; %legend handles

% mark residues that are mutated in comparison with the germline/common ancestor
%===
 ms=6;
 seq0=ma(2,:); % for ch103
 seq0=ma(3,:); % for 3bnc

 if (ii==1)
  seqm=ma(ii,:);
% mind=find(seq0~=seqm);
  mind = find ( (seq0~=seqm).*(seq0~='-').*(seqm~='-') ) ;
  for k=mind
   c=char(clrs(ii)); c=c(1); % color
   plot(k,armsf(k,ii), [c,'o'], 'markersize',ms,'markerfacecolor',c)
   text(k,armsf(k,ii), [seq0(k),seqm(k)])
  end
 end
%===
  end

  xlabel('\it residue');
  ylabel('$\it RMSF(\AA)$', 'interpreter','latex');

% print
 xlim([0 xlims(jj)]);
 legend(hleg,leg);
 set(gcf,'paperpositionmode','auto');
 print(gcf, '-depsc2', ['rmsf',flag,'.eps']);
end

% save data
