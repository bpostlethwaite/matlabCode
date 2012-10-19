function coefs = lsfit (t, b, n)
%
% function coefs = lsfit (t, b, n)
%
% Construct coefficients of the polynomial of
% degree at most n-1 that best fits data (t,b)

t = t(:); b = b(:); % make sure t and b are column vectors
m = length(t);

% long and skinny A
A = ones(m,n);
for j=1:n-1
  A(:,j+1) = A(:,j).*t;
end

% normal equations and solution
B = A'*A; y = A'*b;
coefs = B \ y;