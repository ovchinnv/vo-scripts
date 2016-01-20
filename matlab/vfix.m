% fix diagonals if zero
du=diag(his,1);  du=min(1,du); du=1-du; %ones in place of zeros; elsewhere -- zeros
dl=diag(his,-1); dl=min(1,dl); dl=1-dl;
his2=diag(du,1)+diag(dl,-1) ;
% add matrix with fixes
his=his+his2;
