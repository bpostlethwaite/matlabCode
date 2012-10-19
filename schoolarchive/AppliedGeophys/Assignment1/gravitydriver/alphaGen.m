function [Wmisfit, misfit, alphaG, alphaValue] = alphaGen(G,WtW,D,Mref,n,tol)

% This super long string gets a range value for matrix W, so we can set the
% alpha range correctly (This is because I build them differently than Eldad).
% So I need to switch the ranges automatically depending on the W that is
% fed into the function. 

c =round(abs( full(log10(sum(sum(sum(WtW)))))));
h = waitbar(0,'Optimizing alpha value, please wait');

% Set up sigma two methods
%sigmaD = std(D); % Get standard deviation of data
sigmaD = 0.04*mean(D);
SigN = sigmaD^2*n^2;
% Set up alpha
alpha = logspace(1,7+c,800);
count = 2;
a = round(length(alpha)/2);
MInv = (G*G' + alpha(a).*(WtW))\(G*D + alpha(a).*(WtW)*Mref);
A = (G'*MInv-D)'*(G'*MInv-D);
while  A > SigN + SigN*tol || A < SigN - SigN*tol
    if A > SigN + SigN*tol
        a = round(a - a/(count)) - 1
    else
        a =  round(a + a/(count)) + 1
    end
    count = count +1;
    MInv = (G*G' + alpha(a).*(WtW))\(G*D + alpha(a).*(WtW)*Mref);
    A = (G'*MInv-D)'*(G'*MInv-D); 
    waitbar(count/(0.5*length(alpha)))
end
alphaValue = alpha(a);
close(h)

% Get nice misfit curve
% alphaG = alpha(1:10:length(alpha));
% Wmisfit = zeros(1,length(alphaG));
% misfit = zeros(1,length(alphaG));
% for ii = 1:10:length(alphaG)
%     MInv = (G*G' + alphaG(ii).*(WtW))\(G*D + alphaG(ii).*(WtW)*mref);
%     misfit(ii)=(G'*MInv-D)'*(G'*MInv-D);
%      Wmisfit(ii) = ((W*MInv)'*(W*MInv));
% end
Wmisfit = 0;
misfit = 0;
alphaG = 0;
end