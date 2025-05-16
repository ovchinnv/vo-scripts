% compute & plot BSA
qoct=exist('OCTAVE_VERSION');
if (qoct) ; graphics_toolkit('gnuplot'); end
%
fpath='dcd';
fnames={
'h7-cr8020-noc15_m.sasa'
'h7-cr8020-noc15-fab1.sasa'
'h7-cr8020-noc15-fab2.sasa'
'h7-cr8020-noc15-fab3.sasa'
'h7-cr8020-noc15-notfab1.sasa'
'h7-cr8020-noc15-notfab2.sasa'
'h7-cr8020-noc15-notfab3.sasa'
}

qblock=0; % obtain errors from block averaging or correlation analysis (qblock=0)
% NOTE : qblock=1 gives somewhat higher errors, but both methods seem consistent

%read all data, assuming it fits into the same array (can mod later for generality)
% actually, read the last 10K values (the first 1 or 2 are the minimized pdb coords)
i=1;
numval=7000; % note that this can be set smaller to empirically account for a transient
clear sasa;
for fname=fnames(:)'
 d=load([fpath,'/',fname{1}]);
 sasa(:,i)=d(end-numval+1:end);
 i=i+1;
end

%plot sasa series (they will be noisy and correlated, but we should get a gist of the trends/convergence)
%plot(sasa) ;% all at once ; looks converged, but a very crude vis ;
% use block averaging to compute average + sd
nsasa=numel(fnames);
asasa=mean(sasa,1); % sample means
%sasa=bsxfun(@minus,sasa,asasa); % subtract mean
if (qblock) % for block averaging
for i=1:nsasa
% follow flyvbjerg and petersen 89, p. 463
 d=sasa(:,i);
 j=1;
 bsize(j)=1;
 while 1
  nsamp(j)=numel(d);
%  c0 = sum(d.^2) / nsamp(j); %manual
  c0=var(d,1);
  vsasa(j,i)=c0/(nsamp(j)-1);
  verrsasa(j,i)=vsasa(j,i) / sqrt(2*(nsamp(j)-1));
  if (nsamp(j)<4) ; break ; end % next series will only have one sample, so stop
% otherwise, continue
  j=j+1;
  d = 0.5 * ( d(1:2:end-1) + d(2:2:end) );
  nsamp(j)=numel(d);
  bsize(j)=2*bsize(j-1);
 end
end
% plot estimates of stderr
% note: the error just keeps increasing -- i.e. does not converge ; however, we can still estimate the error
stdsasa=sqrt(vsasa);
%plot(bsize,stdsasa,'o-') ; set(gca, 'xscale','log') ; xlabel('block size') ; ylabel('std err');
%errorbar(sqrt(vsasa),verrsasa,'o-')
%%%%%% index of the blocking data
blkind=10 % by eye
%
stdesasa=stdsasa(blkind,:);
else % qblock
addpath('~/scripts/simonson-util/matlab/');
samps={};
for i=1:nsasa
% follow flyvbjerg and petersen 89, p. 463
 d=sasa(:,i);
 cl=floor(acorr(d)); % correlation length
 samp=[];
 ibeg=1;
 iend=cl;
 while iend<=numel(d)
  samp=[samp mean( d(ibeg:iend) ) ];
  ibeg=iend+1;
  iend=iend+cl;
 end
 samps{i}=samp;
 stdesasa(i)=sqrt(var(samp)/numel(samp));
end
end
%return
% compute & plot BSA :
inds = [1 2 5 ; 1 3 6 ; 1 4 7 ]; % which entries are combines
coefs= 0.5*[-1 1 1 ;-1 1 1 ; -1 1 1 ];
for k=1:size(inds,1)
 ind=inds(k,:);
 coef=coefs(k,:);
 bsa(k) = asasa(ind)*coef';
 bsae(k)= sqrt ( stdesasa(ind).^2 * abs(coef)' );
end

bsa
bsae
data=[bsa;bsae];
save('-ascii', 'h7-cr8020-bsa.dat','data')
% compute average BSA :
