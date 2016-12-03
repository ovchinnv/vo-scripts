% covariance matrices
%
name='3bnc60t';
flag='-ca-hc';
%flag='-ca-lc';

pdbfile=[name,flag,'.pdb'];

struc ;
% define some selections
select ;

% dcd trajectory names :
dcdnames={ [name,flag,'.dcd'] } ;
%
% read dcds :
traj;

% quasiharmonic analysis:
temp=298; %K
[tsquant, tsclass]=calc_quasi(xall, yall, zall, mass, temp)
