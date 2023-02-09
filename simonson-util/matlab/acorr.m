function ac=acorr(cl)
% compute approximate correlation time
  ds=cl;
  dsm=mean(ds);
  dss=std(ds);
  dn=ds-dsm;
  dn=(ds-dsm)/dss;
  dhat=fft(dn);
  dhat2=dhat.*conj(dhat);
  dc=real(ifft(dhat2));
%
%  dc=xcorr(dn,'unbiased'); alternative matlab calc
  dcl=floor(length(ds)/2);
%  dc=autocorr(dn,dcl); % alternative
  dc=dc(1:dcl);
  dc=dc/dc(1) ; % normalization
% find the point at which dc crosses a value near zero
  cutoff= 0.01 ;%empirical
  alpha=1; % empirical integration limit
  k=find(dc<cutoff,1); k=floor(k*alpha) ; 
% integrate the correlation function to this limit, and use the integral to set a correlation length :
  ac=2*trapz(dc(1:k));
end
