% indicate domains on covariance plot

addpath('~/scripts/matlab');
run colours;

ind=[mafwr ; macdr ] ; % defined in malign.mat

legends={ 'FRW1', 'FRW2', 'FRW3', 'FRW4', 'CDR1', 'CDR2', 'CDR3' };

str=ind(:,1); %starting indices
stp=ind(:,2); %stopping indices

l=length(str);

for i=1:l

% draw lines
 off=10;
 x=[1-off, natom+off];
 plot(x,[str(i) str(i)],'color',white,'linewidth',2)
 plot([str(i) str(i)],x,'color',white,'linewidth',2)

% write text labels
 t=text( 0.5 *( str(i) + stp(i) ) , natom+1 , char(legends(i)), 'fontsize',12,'rotation', 90,'fontweight','bold');
% write text labels
 t=text( natom+1, 0.5 *( str(i) + stp(i) ) ,  char(legends(i)), 'fontsize',12,'rotation', 0,'fontweight','bold');
%
end

set(gca, 'position', get(gca,'position')*0.95)