close;

conf='rig';
%conf='pps';

pdb=load(['bfact_',conf,'.dat']);
rms=load(['rmsf_c6',conf,'.dat']);


B=pdb(:,2);
residpdb=pdb(:,1);
rpdb=sqrt(B*3/8)/pi;
% truncate before converter:
%ind=find(residpdb==707);
%residpdb=residpdb(ind:end);
%rpdb=rpdb(ind:end);

rsim=rms(:,3);
residsim=rms(:,2);


%% compute correlation coefficient
%% bring to the same grid;
i=1;
A=[];
ind=0;

for i=1:length(residsim)
 j=find(residpdb==residsim(i));
 if (length(j)>0)
  ind=ind+1;
  A=[ A ; ind rsim(i) rpdb(j)];
 end
end 

%plot(A(:,1), A(:,2),'r');hold on;
%plot(A(:,1), A(:,3),'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=A(:,2);
b=A(:,3);
a=a-mean(a);
b=b-mean(b);

corr=(a'*b)/sqrt((a'*a)*(b'*b))

%%%%%%%%%%%%%% plot %%%%%%%%%%%%
figure;hold on;
%colormap gray
showdom;

plot(residpdb, rpdb, 'k-', residsim, rsim, 'k--', 'Linewidth',1.5);
xlabel('\it Residue ID');
ylabel('\it RMSF (A)');
%legend('calculate\\\\d from experimental B-factors','calculated from MD simulation',2); 
ylim([0 3]);
xlim([707 788]);
box on;

%print(gcf,'-deps',['bfact_c6',conf]);



