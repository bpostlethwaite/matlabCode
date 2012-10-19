% Question 4
clear all
close all
%Hessenberg Matrix

n = 5;
H = triu(rand(n,n),-1);
b = ones(n,1);

%% MAtrix Multiplication way
% tic
% I = eye(n,n);
% M = I;
% MM = I;
% HH = H;
% flip = diag(-1*ones(n-1,1),-1) + eye(n,n);
% for k = 1:n-1
%     M = I;
%     M(k+1,k)  = -HH(k+1,k) / HH(k,k);
%     HH = M*HH;
%     MM = MM * M;
% end
% L = flip.*MM;
% U = HH;
% toc


%% Specific Element Vectorized method
tic
U = H;
for k = 1:n-1
    q = U(k+1,k)/U(k,k);
    U(k+1,k:n) = U(k+1,k:n) - q*U(k,k:n);
    L(k+1,k) = q;
    L(k,k) = 1;
end
L(n,n) = 1;
toc

%% Trials

spy(sparse(L))

y = L\b;
x = U\y;

X = H\b;

X - x;

n = 25;
i = 3;
j = 7;

e = ones(n,1);
S = spdiags([e,e],[-1,0],n,n);

S([i,j],:)=S([j,i],:);
spy(S)
