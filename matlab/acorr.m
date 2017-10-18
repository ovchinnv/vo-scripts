function ac=acorr(cl)
% compute approximate correlation time
  ds=cl;
  dsm=mean(ds);
  dss=std(ds);
  dn=ds-dsm;
  dn=(ds-dsm)/dss;
  dhat=fft(dn);
  dhat2=dhat.*conj(dhat);
  dc=ifft(dhat2);
%
%  dc=xcorr(dn,'unbiased'); alternative matlab calc
  dcl=floor(length(ds)/2);
%  dc=autocorr(dn,dcl); % alternative
  dc=dc(1:dcl);
  dc=dc/dc(1) ; % normalization
% find the point at which dc crosses zero; set block length to alpha times that interval
  alpha = 1.;
  cutoff= 0.0 ;%empirical
  k=find(dc<cutoff,1); k=floor(k*alpha) ; 
  iblock=k;
%  iblock=min(k, floor(length(dn)/2));
  ac=iblock;
end
