
%### converter loops
%#atomselect macro cloop1 { segid CONV and (resid 710 to 711) }
%#atomselect macro cloop2 { segid CONV and (resid 742 to 748) }
%#atomselect macro cloop3 { segid CONV and (resid 759 to 762) }
%#atomselect macro cloop4 { segid CONV and (resid 771 to 773) }
%### converter helices
%#atomselect macro hel1 { segid CONV and (resid 711 to 720) }
%#atomselect macro hel2 { segid CONV and (resid 724 to 729) }
%#atomselect macro hel3 { segid CONV and (resid 730 to 740) }
%#atomselect macro hel4 { segid CONV and (resid 762 to 769) }
%#atomselect macro bets { segid CONV and (( resid 748 to 751 ) or ( resid 754 to 758 ) or ( resid 706 to 710 )) }
%#
%#atomselect macro Cinsert   { segid CONV and (resid 774 to 788) } ;# ok

% places domain separators and domain labels;
ind=[ 742, 748, ... % cl2
759, 762, ... % cl3
770, 773, ... % cl4
711, 720, ... % hel1
724, 729, ... % hel2
730, 740, ... % hel3
762, 769, ... % hel4
774, 788 ... % CI
];

% custom colors
colors=[ [255 255 255]*0.8;... %lp2
[255 255 255]*0.8;... %lp3
[255 255 255]*0.8;... % lp4
235 0 55;... %H1
0 212 9;... %H2
98 98 147;... %H3
243 0 160;... %H4
185 111 111]; %cinsert
colors=colors/256;

%colors=colors.^0.7; % fade

l=length(ind);
ind2=reshape(ind,2, l/2)';

str=ind2(:,1); %starting indices
stp=ind2(:,2); %stopping indices

%plot the lines

l=length(str);
yl=10;

c=colormap;
cinc=length(c)/l;

gray=[0.8 0.8 0.8];
white=[1 1 1];
red=[1 0 0];

for i=1:l
 x=linspace(str(i), stp(i),2);
% h=rectangle('position', [ x(1) 0 x(2)-x(1) ylim ], 'facecolor', gray, 'edgecolor',[0 0 0]);
% h=rectangle('position', [ x(1) 0 x(2)-x(1) yl ], 'facecolor', c(cinc*i,:), 'edgecolor',[0 0 0]);

% prescribed colors:
 h=rectangle('position', [ x(1) 0 x(2)-x(1) yl ], 'facecolor', colors(i,:), 'edgecolor',[0 0 0]);

end  
%box on;


% print;

%print(gcf, '-deps2', 'domains.ps');
%print(gcf, '-painters', '-dpsc2', 'domains.ps');

