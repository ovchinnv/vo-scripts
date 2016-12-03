
% places domain separators and domain labels;
%sw1=[36, 47];
%sw2=[66, 85];
%ploop=[18, 25];
%intersw=[48, 65];
%loop1=[111, 119] ; % flexible loop that interacts with sw 2
%
%ind=[sw1 ; sw2 ; ploop ; intersw ;loop1 ];

ind=[mafwr ; macdr ] ; % assume that we have read malign.mat

legends={ 'FRW1', 'FRW2', 'FRW3', 'FRW4', 'CDR1', 'CDR2', 'CDR3' };

addpath('~/scripts/matlab');
run colours
colors=[pink ; pink; pink ; pink ; steelblue ; steelblue ; steelblue ; steelblue ];

%colors=colors.^0.7; % fade

str=ind(:,1); %starting indices
stp=ind(:,2); %stopping indices

%plot the lines

l=length(str);
yl=6;

cols=colormap;
cinc=length(cols)/l;

gray=[0.8 0.8 0.8];
white=[1 1 1];
red=[1 0 0];

for i=1:l
 x=linspace(str(i), stp(i),2);
% h=rectangle('position', [ x(1) 0 x(2)-x(1) ylim ], 'facecolor', gray, 'edgecolor',[0 0 0]);
% h=rectangle('position', [ x(1) 0 x(2)-x(1) yl ], 'facecolor', c(cinc*i,:), 'edgecolor',[0 0 0]);

% prescribed colors:
 h=rectangle('position', [ x(1) 0 x(2)-x(1) yl ], 'facecolor', colors(i,:), 'edgecolor',[0 0 0]);
% write text labels
 t=text( 0.5 * ( x(2) + x(1) ) + 0 , 3.25 , char(legends(i)), 'fontsize',12,'rotation', 90,'fontweight','bold');
%
end
%box on;
% print;
%print(gcf, '-deps2', 'domains.ps');
%print(gcf, '-painters', '-dpsc2', 'domains.ps');

