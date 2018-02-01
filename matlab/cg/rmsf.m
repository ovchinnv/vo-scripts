%

fname_c='rmsf_clpx_cgs.dat';
fname_f='rmsf_cgsmall.dat';

rf=load(fname_f);
rc=load(fname_c);

figure; hold on;

i=rf(:,1);
ires=rf(:,2);
rms=rf(:,3);

plot(rf(:,1),rf(:,3),'k.-');
plot(rc(:,1),rc(:,3),'r.-');

%plot(ires,rms,'k*-');
%bar(ires,rms,'k*-');




