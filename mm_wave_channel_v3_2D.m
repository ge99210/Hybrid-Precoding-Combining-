function [H_n, Phi_AOD, Phi_AOA, Alpha] = mm_wave_channel_v3_2D(Nt, Nr, Nc, Np, sig)
% mmWave channel model with Np rays per cluster for ULA in the y-axis
% Georgios K. Papageorgiou 19/07/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nt: number of transmit antennas
% Nr: number of receive antennas
% Nc: number of clusters 
% Np: number of paths (rays) per cluster
% sig: the standard deviation of the random angle from the center of the cluster 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H_n : the channel (normalized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the steering vector function of a ULA of lambda/2 spacing along the y-axis (considering that theta = pi/2 
% in ISO Spherical coordinate system)
a = @(phi,N) exp(-1j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);
%
% The azimuth angles for the cluster's centers
Phi_AOD_m = -pi+2*pi*rand(Nc,1);
Phi_AOA_m = -pi+2*pi*rand(Nc,1);
Alpha = (1/sqrt(2))*(randn(Nc,Np)+1j*randn(Nc,Np)); % CN(0,1) complex Gaussian with unit variance 

% Initialization
Phi_AOD = zeros(Nc,Np);
Phi_AOA = zeros(Nc,Np);
H = zeros(Nr,Nt);
 for i=1:Nc
     % The final azimuth angles from each cluster i    
     phi_AOD = tr_laprnd(Np, 1, Phi_AOD_m(i), sig);
     phi_AOA = tr_laprnd(Np, 1, Phi_AOA_m(i), sig);    
         for l=1:Np
         Phi_AOD(i,l) =  phi_AOD(l);
         Phi_AOA(i,l) =  phi_AOA(l);
         alpha = Alpha(i,l);
         A_t = a(phi_AOD, Nt); % Tx steering vector in 2D
         A_r = a(phi_AOA, Nr); % Rx steering vector in 2D

         H = H + alpha*A_r*A_t'; % channel
         end
 end
 H_n = sqrt(Nt*Nr/(Nc*Np))*H;
end

function y  = tr_laprnd(m, n, mu, sigma)
%Truncated LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution
%   Author  : 
%   Date    : 
%Check inputs
if nargin < 2
    error('At least two inputs are required');
end
if nargin == 2
    mu = 0; sigma = 1;
end
if nargin == 3
    sigma = 1;
end
% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b*sign(u).*log(1- 2*(1-exp(-pi/b))*abs(u));
end
