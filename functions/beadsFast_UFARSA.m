function [x, f] = beadsFast_UFARSA(y, fc, lam0)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code provides a fast version of the BEADS method, by adapting its original 
% MATLAB implementation written by Laurent Duval (doi: 10.1016/j.chemolab.2014.09.014):
% https://de.mathworks.com/matlabcentral/fileexchange/49974-beads--baseline-estimation-and-denoising-w--sparsity--chromatogram-signals-
% 
% IMPORTANT: We don't claim that the results of this fast-implementation 
% are exactly the same as the original code. However, we found that this 
% fast version can generate qualitatively similar results, which were 
% sufficient to remove slowly varying drifts from our calcium imaging data.
% Therefore, to use BEADS method out of UFARSA, please use the orignal 
% code (see the link above).
% In addition to our fast-version implementation, we made several assumptions about the parameters: 
% lam1=lam2=1; d =1; pen = 'L1_v2'; r=1; Nit=1; EPS0 = 0.6; EPS1 = 1e-2; no "cost" output.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R[x, f] = beadsFast_UFARSA(y, fc, lam0)
%
% INPUT
%   y: Noisy observation
%   fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%   lam0: Regularization parameter
%
% OUTPUT
%   x: Estimated sparse-derivative signal
%   f: Estimated baseline (i.e. the trace of slowly varying drifts)
%
% Author: Vahid Rahmati (December, 2017)

%%
EPS0 = 0.6;    % cost smoothing parameter for x 
EPS1 = 1e-2;   % cost smoothing parameter for derivatives 
r    = 1;      % Asymmetry ratio

y = y(:);
x = y;
N = length(y);

[A, B] = BAfilt_fast(fc, N);

e = ones(N-1, 1);
D1 = spdiags([-e e], [0 1], N-1, N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D = [D1;  D2];
BTB = B'*B;

w = ones(2*N-3,1);
d = BTB * (A\y);

gamma = ones(N, 1);
Lambda = sparse(1:2*N-3, 1:2*N-3, w.*(1./( abs(D*x) + EPS1)), 2*N-3, 2*N-3);

k = abs(x) > EPS0;
gamma(~k) = ((1 + r)/4) / abs(EPS0);
gamma(k) = ((1 + r)/4) ./  abs(x(k));
Gamma = sparse(1:N, 1:N,gamma,N,N);

M = 2 * lam0 * Gamma + D' * Lambda * D;
x = A * ((BTB + A'*M*A)\d);

f = y - x - B*(A\(y-x));

end

%% local function
function [A, B] = BAfilt_fast(fc, N)
% [A, B] = BAfilt(d, fc, N)
%
% Computing banded matrices for zero-phase high-pass filter.

d = 1;
b = [-1 2 -1];
omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = [1  2 1];
a = b + t*a;
   
A1 = sparse(1:N,1:N,a(2),N,N);
A2 = sparse(2:N,1:N-1,a(1),N,N);
A  = A2 + A1 + A2'; % A: Symmetric banded matrix
 
B1 = sparse(1:N,1:N,b(2),N,N);
B2 = sparse(2:N,1:N-1,b(1),N,N);
B  = B2 + B1 + B2'; % B: banded matrix

end

