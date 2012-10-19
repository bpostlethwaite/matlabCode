function tab = romberg (r,a,b)
%
% Generate the Romberg table
%
% The integrand is defined by the following,
% even though this is not cool

  h0 = b-a;
  h(1) = h0;
  for j=2:r
    h(j) = h(j-1)/2;
  end
  tab = zeros(r,r);

  for j=1:r
    x = a:h(j):b;
    ex = exp(-x.^2);
    tab(j,1) = h(j)*sum(ex);
    tab(j,1) = tab(j,1) - h(j)/2*(ex(1)+ex(length(x)));
  end

% Higher order approximations

  for j=2:r
    for k=2:j
      tab(j,k)=tab(j,k-1)+ (tab(j,k-1)-tab(j-1,k-1))/(4^(k-1)-1);
    end
  end