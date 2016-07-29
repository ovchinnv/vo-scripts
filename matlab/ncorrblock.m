%% ncorrsamp
%% use correlation to approximate the number of independent samples in correlated data
%%
function nsamp=ncorrblock(x)
  alpha=1;
  dn=x-mean(x);
  dhat=fft(dn);
  dhat2=dhat.*conj(dhat);
  dc=ifft(dhat2);
  dc=dc(1:floor(length(dc)/2)) ;
  dc=dc/dc(1) ; % normalize
% find the point at which dc crosses zero; set block length to alpha times that interval
  k=find(dc<0,1); k=floor(k*alpha) ; iblock=min(k, floor(length(dn)/2));
  nsamp=length(x)/iblock ; %effective number of independent samples
end
