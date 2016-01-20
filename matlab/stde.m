%
% computation of standard error using method in Rapaport, p. 86
% note that this estimate agrees with what I get out of Kunsch 98
%
% sample size

x=pbe1(:);

n=length(x);

% arrays containing block averages
xb=zeros(n,1);
% arrays containing block variances
xbv=zeros(n,1);

% variance estimates 
blocks=ceil(log2(n));
varx=zeros(blocks,1);
varxv=zeros(blocks,1);
bsize=zeros(blocks,1)

l=0;
d=1; % block size
while (2*d<=n)
d=2^l;
l=l+1;
% compute block averages:
fprintf(['Block size: ',num2str(d),'\n']);

k=0;
for i=1:d:n
 k=k+1;
 j=min(i+d-1,n);
 xb (k) = ( mean(x(i:j)) )     ; % block-average 
 xbv(k) = ( var (x(i:j))  )    ; % block-variance
end

%varp(l) = sum ( peb(1:k).^2 - pea^2 ) / ( k )
%varp(l) = sum ( peb(1:k).^2 - pea^2 ) / ( k-1 )
varx(l)  = var(xb(1:k));

varxv(l) = var(xbv(1:k));

bsize(l) = d;

end

semilogx(bsize(1:l),sqrt(varx(1:l)));
xlabel('Block size');

vfig=figure(2)
semilogx(bsize(1:l),sqrt(varxv(1:l)));
xlabel('Block size');


[sqrt(varx) , sqrt(varxv)]

