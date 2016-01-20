%
close all;

if (~exist('read'))
 read=1;
end

if (read==1)
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'curv5.dat' 'curv6.dat' 'curv7.dat' 'curv8.dat' 'curv9.dat' 'curv10.dat' 'curv11.dat' 'curv12.dat' 'curv13.dat' ...
 'curv14.dat' 'curv15.dat' 'curv16.dat' 'curv17.dat' 'curv18.dat'  'curv19.dat' 'curv20.dat' 'curv21.dat' 'curv22.dat' ...
 'curv23.dat' 'curv24.dat' 'curv25.dat' 'curv26.dat' 'curv27.dat'  'curv28.dat' 'curv29.dat' 'curv30.dat' 'curv31.dat' 'curv32.dat' ...
 'curv33.dat' 'curv34.dat' 'curv35.dat' 'curv36.dat' 'curv37.dat' 'curv38.dat' ...
 'curv39.dat' 'curv40.dat' 'curv41.dat' 'curv42.dat' 'curv43.dat' 'curv44.dat' ...
 'curv45.dat' 'curv46.dat' 'curv47.dat' 'curv48.dat' 'curv49.dat' 'curv50.dat' 'curv51.dat' 'curv52.dat' ...
 'curv53.dat' 'curv54.dat' 'curv55.dat' 'curv56.dat' 'curv57.dat' 'curv58.dat' 'curv59.dat' 'curv60.dat' ...
 'curv61.dat' 'curv62.dat' 'curv63.dat' 'curv64.dat' 'curv65.dat' ...
 'curv66.dat' 'curv67.dat' 'curv68.dat' 'curv69.dat' 'curv70.dat' ...
 'curv71.dat' 'curv72.dat' 'curv73.dat' 'curv74.dat' 'curv75.dat' 'curv76.dat' 'curv77.dat' 'curv78.dat' 'curv79.dat' ...
 'curv80.dat' 'curv81.dat' ...
 }
%
clear crv;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  crv=load(fname);
 else
  crv=[crv; load(fname)];
 end 
end
read=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[niter,nstring]=size(crv);

ind=[1:niter];
d=20;
crvsum = sum(crv(:,2:end),2)
crvsum = smooth2(ind, crvsum, d);
%
%
plot(ind, crvsum );

box on;
ylabel('Replica index');
xlabel('iteration');

set(gcf, 'paperpositionmode', 'auto');

%print(gcf, '-depsc2', 'slen.eps');
%print(gcf, '-djpeg100', 'slen.jpg');



