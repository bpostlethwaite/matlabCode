clear all; close all
% Play with matrices

%{
G = randn(5,5);
B = orth(G);
subplot(1,2,1)
imagesc(G)
subplot(1,2,2)
imagesc(B)
%}

S = zeros(20,60)

for ii = [10,20,40,50]
S(:,[1,2,3]+ii) = 1;

end

spy(S)
title('Block-Sparsity with d = 3')