% Dave Mihn's version
% ----
% WHAM
% ----
%
% Calculates the from a set of biased simulations using the 
% weighted histogram analysis method, as described by Roux.
% 
% B. Roux. The calculation of the potential of mean force using computer
% simulations. Computer Physics Communications 91, 275 (1995). 
%
% [rU,xspace] = WHAM(Xs,Bxs,a,xspace)
%
% Parameters
% Xs:      a WxS matrix containing X coordinate values,
%          where W is the window and S is the step
% Bxs:     a Wx1 vector containing the location of the 
%          umbrella window centers.
% xspace:  a 1xB vector marking the edges of the bins along X,
%          where B is the number of bins
%
% Outputs
% rU:      the PMF
% xspace:  a 1xB vector marking the edges of the bins along X,
%          where B is the number of bins
%
function [rU,xspace] = WHAM(Xs,Bxs,xspace)

KbT = 1.3806503E-2*300;     %%% energy unit, (pN nm)
ks = 0.1;                   %%% bias spring constant, pN/nm

%%% Default range for X
if nargin < 3
    xspace = linspace(min(Xs(1:end)),max(Xs(1:end)),50);
end
xsp = xspace(2)-xspace(1);


steps = size(Xs(1:end),2);  %%% total number of data points
numsims = size(Bxs,1);      %%% number of windows
numbins = size(xspace,2);   %%% number of bins

%%% Histograms of data
Hx = histc(Xs',xspace)/steps;

%%% WHAM variables
Fx_old = zeros(1,numsims);   %%% Free energy constants Fx_i   
Px = zeros(1,numbins);       %%% Probabilities from equation (8)

%%% Variables tracking the progress of change
Fprog = [];
iter = 0;
change = 0.05;

%%% WHAM iterations
while change > 0.001
    numPx = zeros(1,numbins); %%% Numerator for equation (8)
    denPx = zeros(1,numbins); %%% Denominator for equation (8)
    Px = zeros(1,numbins);    %%% Probabilities from equation (8)
    Fx = zeros(1,numsims);    %%% Free energy constants F_i

    %%% Calculate probabilities from free energies
    for sim = 1:numsims
        numPx = numPx + Hx(:,sim)';
        Ubias = 0.5*ks*((Bxs(sim)-(xspace+xsp/2)).^2);
        denPx = denPx + sum(Hx(:,sim))*exp((Fx_old(sim)-Ubias)/KbT);
    end
    Px = numPx./denPx;
    for sim = 1:numsims
        Ubias = 0.5*ks*((Bxs(sim)-(xspace+xsp/2)).^2);
        Fx(sim) = sum(Px.*exp(-Ubias/KbT));
    end
    Fx = -KbT*log(Fx);
    Fx(2:end) = Fx(2:end)-Fx(1);

    Fx_old = Fx;
    
    %%% Measure change in free energy constants
    Fprog = [Fprog; Fx];
    if iter > 2
        change = max(abs((Fprog(end-1,2:end)-Fprog(end,2:end))));
    end
    iter = iter + 1;
end

%%% Set minima to zero
rU = -log(Px); 
rU = rU - min(rU(1:end));
