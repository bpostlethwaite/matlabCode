function demo(interactive)
% DEMO Demonstrates the use of the SPGL1 solver.
%
% See also SPGL1.
    
%   demo.m
%   $Id: demo.m 315 2007-08-03 18:17:18Z mpf $
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
    if nargin < 1 || isempty(interactive), interactive = true; end
    
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);
    
    % Create random m-by-n encoding matrix and sparse vector
    m = 50; n = 128; k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n);
    x0 = zeros(n,1); x0(p(1:k)) = randn(k,1);

    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
    [x,r,g,info] = spgl1(A, b, tau, [], []);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\nPress <return> to continue ... \n');

    
    if interactive, pause; end
    
    
    % -----------------------------------------------------------
    % Solve the basis pursuit (BP) problem:
    %
    %    minimize ||x|_1 subject to Ax = b
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit (BP) problem:\n');
    fprintf('%%                                      \n');
    fprintf('%%   minimize ||x||_1 subject to Ax = b \n');
    fprintf('%%                                      \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    % Set up vector b, and run solver 
    b     = A * x0;  % Signal
    tau   = 0;       % Initial value of tau
    sigma = 0;       % Desired ||Ax - b||_2
    opts  = spgSetParms('iterations',250,'verbosity',1);
    [x,r,g,info] = spgl1(A, b, tau, sigma, [], opts);
    
    figure(1); subplot(2,2,1);
    plot(1:n,x,'b', 1:n,x0,'ro');
    legend('Recovered coefficients','Original coefficients');
    title('(a) Basis Pursuit');

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(a).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\nPress <return> to continue ... \n');

    
    if interactive, pause; end

    
    % -----------------------------------------------------------
    % Solve the basis pursuit denoise (BPDN) problem:
    %
    %    minimize ||x|_1 subject to ||Ax - b||_2 <= 0.1
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit denoise (BPDN) problem:  \n');
    fprintf('%%                                                  \n');
    fprintf('%%   minimize ||x||_1 subject to ||Ax - b||_2 <= 0.1\n');
    fprintf('%%                                                  \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    % Set up vector b, and run solver
    b = A * x0 + randn(m,1) * 0.075;
    tau   = 0;          % Initial value of tau
    sigma = 0.10;       % Desired ||Ax - b||_2
    xs    = zeros(n,1); % Initial guess for x
    opts = spgSetParms('iterations',1000,'optTol',1e-4,'verbosity',1);
    [x,r,g,info] = spgl1(A, b, tau, sigma, xs, opts);
    
    figure(1); subplot(2,2,2);
    plot(1:n,x,'b', 1:n,x0,'ro');
    legend('Recovered coefficients','Original coefficients');
    title('(b) Basis Pursuit Denoise');

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(b).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\nPress <return> to continue ... \n');


    if interactive, pause; end


    % -----------------------------------------------------------
    % Solve the basis pursuit (BP) problem in COMPLEX variables:
    %
    %    minimize ||x|_1 subject to Ax = b
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the basis pursuit (BP) problem in COMPLEX variables:\n');
    fprintf('%%                                                    \n');
    fprintf('%%   minimize ||x||_1 subject to Ax = b               \n');
    fprintf('%%                                                    \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    % Create partial Fourier operator with rows idx
    idx = randperm(n); idx = idx(1:m);
    opA = @(x,mode) partialFourier(idx,n,x,mode);

    % Create sparse coefficients and b = 'A' * x_0;
    x0c = zeros(n,1);
    x0c(p(1:k)) = randn(k,1) + sqrt(-1) * randn(k,1);
    b = opA(x0c,1);
    
    tau   = 0;
    sigma = 0;
    opts = spgSetParms('iterations',1000,'optTol',1e-4,'verbosity',1);
    [x,r,g,info] = spgl1(opA,b,tau,sigma,[],opts);
    
    figure(1); subplot(2,2,3);
    plot(1:n,real(x),'b+',1:n,real(x0c),'bo', ...
         1:n,imag(x),'r+',1:n,imag(x0c),'ro');
    legend('Recovered (real)', 'Original (real)', ...
           'Recovered (imag)', 'Original (imag)');
    title('(c) Complex Basis Pursuit');
    
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(c).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\nPress <return> to continue ... \n');

    
    if interactive, pause; end
    
    
    % -----------------------------------------------------------
    % Sample the Pareto frontier at 100 points:
    %
    %    phi(tau) = minimize ||Ax-b||_2 subject to ||x|| <= tau
    %
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Sample the Pareto frontier at 100 points:\n');
    fprintf('%%                                              \n');
    fprintf('%%   phi(tau) = minimize ||Ax-b||_2 subject to ||x|| <= tau\n');
    fprintf('%%                                              \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('\nComputing sample');
    
    % Set up vector b, and run solver    
    b   = A*x0;
    x   = zeros(n,1);
    tau = linspace(0,1.05 * norm(x0,1),100);
    phi = zeros(size(tau));
    
    opts = spgSetParms('iterations',1000,'verbosity',0);
    for i=1:length(tau)
        [x,r,g,info] = spgl1(A,b,tau(i),[],x,opts);
        phi(i) = info.rNorm;
        if ~mod(i,10), fprintf('...%i',i); end
    end
    fprintf('\n');
    
    figure(1); subplot(2,2,4);
    plot(tau,phi);
    title('(d) Pareto frontier');
    xlabel('||x||_1'); ylabel('||Ax-b||_2');
    
    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(d).\n');
    fprintf([repmat('-',1,80), '\n']);
    
end % function demo    

function y = partialFourier(idx,n,x,mode)
    if mode==1
       % y = P(idx) * FFT(x)
       z = fft(x) / sqrt(n);
       y = z(idx);
    else
       z = zeros(n,1);
       z(idx) = x;
       y = ifft(z) * sqrt(n);
    end
 end % function partialFourier
