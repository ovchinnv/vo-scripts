% octave file to load acemd output
% to use Matlab need to throw away lines beginning with # which octave treats as comments
% this can be done with grep, e.g. grep -ve "^[TCL|#]" [oldfile] > [newfile]
% which will get rid of lines starting with # or with TCL
%

f='npt5eq'; % this is the log

cmd=[" grep -ve '^[TCL|#]' ", f," > temp_"]

system ( cmd )

%d=load(f);
d=load('temp_');

% these are the components from the log file :
%     Step      Bond     Angle     Dihed      Elec       VDW        PE        KE  Externa     Total      Temp      Pres   PresAve
%
step=1;
bond=2;
angl=3;
dihe=4;
elec=5;
vdw=6;
pe=7;
ke=8;
exte=9;
total=10;
temp=11;
pres=12;
prav=13;

dt=2 ; % time step in fs
time=d(:,step) * 2 / 1000000 ; % time in ns 


plot(time, d(:,pres),'k--') ; hold on
plot(time, d(:,prav),'r-');

legend('instantaneous','moving average');

xlabel('time(ns)');
ylabel('Pressure (bar) ');

ylim([-10,10]); % set axis limits

% average pressure : 
% discard initial portion of trajectory
n=length(time) ;
plongave=mean(d( floor(n*0.1) : end, prav)); %
fprintf(['Average pressure : ',num2str(plongave),'\n'])

