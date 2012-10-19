function [x, itn] = oneProjectorSort(b, tau)
% [x, itn] = oneProjectorSort(b, tau) 
% Return the orthogonal projection of the vector b onto the L1 ball.
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||x||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%


% Initialization
n     = length(b);
x     = zeros(n,1);
bNorm = norm(b,1);

% Check for quick exit.
if (tau >= bNorm), x = b; itn = 0; return; end
if (tau <  eps  ),        itn = 0; return; end

% Preprocessing
s         = sign(b);
b         = abs(b);
[b,idx]   = sort(b,'descend'); % Descending.

csb       = -tau;
alphaPrev = 0;
 for j= 1:n
    csb       = csb + b(j);
    alpha     = csb / j;
    
    % We are done as soon as the constraint can be satisfied
    % without exceeding the current minimum value of b
    if alpha >= b(j)
       break;
    end
    
    alphaPrev = alpha;
 end

% Set the solution by applying soft-thresholding with
% the previous value of alpha
x = max(0,b - alphaPrev);

% Postprocessing
x(idx) = s(idx) .* x;
itn    = j;
