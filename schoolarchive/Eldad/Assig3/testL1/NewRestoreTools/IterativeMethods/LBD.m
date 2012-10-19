function [U, B, V] = LBD(A, U, B, V, P, options)
%
%     [U, B, V] = LBD(A, U, B, V, P, options)
%
%  Perform one step of Lanczos bidiagonalization with or without
%  reorthogonalization, no preconditioner here.
%
% Input:
%        A - matrix
%     U, V - accumulation of vectors
%        B - bidiagonal matrix
%        P - input ignored
%  options - structure from HyBR (see HyBRset)
% Output:
%     U, V - updated "orthogonal" matrix
%        B - updated bidiagonal matrix
%
%  Refs: 
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%
reorth = strcmp(HyBRget(options,'Reorth'), {'on'});

[m, k] = size(U);

if k == 1
  v = A'*U(:,k);
else
  v = A'*U(:,k) - B(k, k-1)*V(:,k-1);
  
  if reorth %  the next step is MGS reorthogonalization
    for j = 1:k-1
      v = v - (V(:,j)'*v)*V(:,j);
    end
  end
  
end
alpha = norm(v);
v = v / alpha;

u = A*v - alpha*U(:,k);

  if reorth %  the next step is MGS reorthogonalization
    for j = 1:k
      u = u - (U(:,j)'*u)*U(:,j);
    end
  end
  
beta = norm(u);
u = u / beta;
  
U = [U, u];
V = [V, v];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];

