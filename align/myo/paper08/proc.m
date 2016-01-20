
nseq=38;

[id, seq, nres, tag]=textread('align2.txt', '%s %s %d %s','delimiter',' ');

nlines=length(id);
neol=nlines/nseq;

id=id(1:nseq);
tag=tag(1:nseq);

for i=1:nseq
 for j=1:neol-1
  seq(i)={ [char(seq(i)), char(seq(i+nseq*j))] };
 end
end
seq=seq(1:nseq);


% convert sequence into character array
seqc=char(seq);

% generate residue numbers for each sequence

resmax=length(char(seq(1)));
resid=zeros(nseq, resmax);

%fprintf('%s\n\n', 'Numbering sequences...');
for i=1:nseq
% fprintf('%d\n',i);
 for j=1:resmax
%
   if (seqc(i,j) ~= '-' )
    k=j-1;
    resid(i,j)=1;
    while (k>1) 
     if (resid(i,k)>0)
      resid(i,j)=resid(i,k)+1;
      break;
     end
     k=k-1;
    end
   else
    resid(i,j)=-1;
   end
%
 end
end


%%%%%%%%%%% now, can produce graphical alignment %%%%%%%%%%%%%%
% 1) put sequences in desired order
order=[20 21 22 17 4 12 25 10 35 26 27 28 30 32 33 34 13 14 15 5 6 7 8 23 24 11 16 1 18 19 37 38 9 2 3 31 36 29];
%reverse lookup
order2=zeros(size(order));
for i=1:nseq
 ind=find(order==i);
 order2(i)=ind(1);
end

%tag(order2) prints tags in order from myo1 to myo18


%%%%%% specify MV residues and print equivalent residues in other seqs %%%%
% also indicate how close a match is achieved


pos=['RK']; % 2
hphob=['AILMFWYV'];%10
pol=['NCQGHPST']; %18
neg=['DE']; %20


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% residues in MV sequence
res=[246 639 398 566 405 570 ... % cleft OP
...%     174 204 178 193 ... % HF/HH/HG
...%     195 ... % HG/HH
     197 202 203 206 412 423 ... %nonpolar interface between HO/HH/HG
     219 442 ... %sw1/2 salt bridge
     449 677 680 684 ] %relay/HW connunication
      
%%%%% note row 17 corresponds to chicken MV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% produce alignment graph %


n=length(res);
nalign=length(order2);
residues=char(zeros(n,nalign));

resnum=zeros(n,nalign);

for i=1:n
% find index in sequence
 ind=find(resid(17,:)==res(i));
 resnum(i,:)=resid(order2,ind)';
 residues(i,:)=seqc(order2,ind)';
end
tags=tag(order2);
ids=id(order2);

%residues
% plot !
close;
f=figure('position',[10, 20, 900 30*nalign]); hold on;

dy=-15; 
dx=40;
yl=abs(dy) * ( nalign + 1);
xl=900;
axis([0 xl 0 yl]); box on;

x0=180; y0=abs(nalign*dy);


% first, print tags
%
for j=1: nalign
  text(10,y0+(j-1)*dy,[char(tags(j)),'(',char(ids(j)),')'], 'FontName','courier','FontSize',8,'FontWeight','normal', 'interpreter','none');
end

for i=1: n
 m5res=residues(i,order(17));
 for j=1: nalign
  if (order2(j)==17) 
   c='black';w='bold';
  elseif (m5res==residues(i,j))
   c='black';w='normal';
  elseif ( ( (sum (pos==m5res)==1) && (sum(pos==residues(i,j))==1)) || ...
        ( (sum (neg==m5res)==1) && (sum(neg==residues(i,j))==1)) || ...
        ( (sum (pol==m5res)==1) && (sum(pol==residues(i,j))==1)) || ...
        ( (sum (hphob==m5res)==1) && (sum(hphob==residues(i,j))==1)) )
     c='blue';w='demi';
  else
   c='red';w='demi';
  end
% length of number field
  resno=num2str(resnum(i,j));
  while ( length(resno) < 4 )
   resno=[' ', resno];
  end 
%
  text(x0 + (i-1) * dx, y0 + (j-1) * dy, [resno, residues(i,j)], 'FontName','courier','FontSize',8,'FontWeight',w,'color',c);
 end
end

set(gca,'Ytick', []); set(gca,'YtickLabel',[] );
set(gca,'Xtick', []); set(gca,'XtickLabel',[] );

% print
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-depsc', 'align.eps');
print(gcf, '-djpeg100', 'align.jpg');

