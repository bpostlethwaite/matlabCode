function [npar nperp ts tn] = degreeshift2(tor)

degree = 0:1:89;
for ii = 1:length(degree);
npar(:,ii) = [tand(degree(ii));1]./sqrt( tand(degree(ii))^2 +1);
nperp(:,ii) = [npar(2,ii);-npar(1,ii)];
t(:,ii) = tor*nperp(:,ii);
ts(ii)  = t(:,ii)'*npar(:,ii);
tn(ii) =  t(:,ii)'*nperp(:,ii);
end

npar(:,91) = [1;0];
nperp(:,91)=[0;-1];
t(:,91) = tor*nperp(:,91);
ts(91)  = t(:,91)'*npar(:,91);
tn(91) =  t(:,91)'*nperp(:,91);

for ii = length(degree) +2 : 171
npar(:,ii) = npar(:,182-ii).*[1;-1];
nperp(:,ii) = [npar(2,ii);-npar(1,ii)];
t(:,ii) = tor*nperp(:,ii);
ts(ii)  = t(:,ii)'*npar(:,ii);
tn(ii) =  t(:,ii)'*nperp(:,ii);
end




end