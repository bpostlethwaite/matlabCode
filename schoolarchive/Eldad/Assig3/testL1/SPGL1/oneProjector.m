function [x, itn] = oneProjector(b, tau)
% [x, itn] = oneProjector(b, tau) 
% Return the orthogonal projection of the vector b onto the L1 ball.
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||x||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1.

%   oneProjector.m
%   $Id: oneProjector.m 272 2007-07-19 05:40:29Z mpf $
%
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------

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
