
nmod=200;
ncol=10; % 2 4 6 8 10
nrow=40; % 1 6 11 ...

proj=zeros(nmod,nmod);

flag='mde_mass';

d = load(['dot',flag,'.dat']);

for j=1:200
 for i=1:5
%  proj(1,i:5:1+195)=d(1:40,2*i);
  proj(j,i:5:i+195)=d(40*(j-1)+1:40*(j-1)+40,2*i)';
 end
end

% 1:100 -- rigor modes
% 101:200 -- postrigor modes

%proj1=proj(100:-1:1,200:-1:101);
proj1=proj(100:-1:91,200:-1:191);

%pcolor(proj1); caxis([-1 1]);
pcolor(abs(proj1)); caxis([0 1]);

shading flat;colorbar;
xlabel('\it rigor modes');
ylabel('\it postrigor modes');
title('Scalar product of rigor/postrigor PCA modes (lowest 100 modes)');

box on;
print(gcf,'-depsc',['pca_dot_',flag,'10.eps']);
