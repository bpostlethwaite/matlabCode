function [Q,T] = lanczos(A,Q,k)
%
% function [Q,T] = lanczos(A,Q,k)
%

% preallocate for speed
alpha=zeros(k,1); beta=zeros(k,1);

Q(:,1) = Q(:,1)/norm(Q(:,1)); 
beta(1,1)=0;

for j=1:k
  w=A*Q(:,j);
  if j>1
    w=A*Q(:,j)-beta(j,1)*Q(:,j-1);
  end
  alpha(j,1)=Q(:,j)'*w;
  w=w-alpha(j,1)*Q(:,j);
  beta(j+1,1)=norm(w);

  if abs(beta(j+1,1))<1e-10 
    disp('Zero beta --- returning.');
    T=spdiags([beta(2:j+1) alpha(1:j) beta(1:j)],-1:1,j+1,j);
    return
  end
  Q(:,j+1)=w/beta(j+1,1);
end
T=spdiags([beta(2:end) alpha beta(1:end-1)],-1:1,k+1,k);