% indicate domains on covariance plot

addpath('~/scripts/matlab');
run colours;

ind=[mafwr ; macdr ] ; % defined in malign.mat

legends={ 'FRW1', 'FRW2', 'FRW3', 'FRW4', 'CDR1', 'CDR2', 'CDR3' };

str=ind(:,1); %starting indices
stp=ind(:,2); %stopping indices

% note that the above domain definitions are with respect to 3bnc sequence
% need to convert to unwrapped/aligned sequence :
seqind=3 ;% 3bnc60 gl
seqthis=find(ismember(files,name(1:end-1))); % end-1 to drop 't' et the end
inda=zeros(size(ind)) ;
for i=1:length(ind)
 for j=1:2
  k=find(maind(:,seqind)==ind(i,j));
  inda(i,j)=maind(k,seqthis);
 end
end
%
str=inda(:,1); %starting indices
stp=inda(:,2); %stopping indices

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
