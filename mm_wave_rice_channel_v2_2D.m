function [H_n, Phi_AOD, Phi_AOA, Alpha] =...
    mm_wave_rice_channel_v2_2D(Nt, Nr, Nc, Np, sig, k)
% mmWave Rician channel model with Np rays per cluster for ULA in the
% y-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required MATLAB file: mm_wave_channel_v1_2D.m
% Nt: number of transmit antennas
% Nr: number of receive antennas
% Nc: number of clusters 
% Np: number of paths (rays) per cluster
% sig: the standard deviation of the random angle from the center of the cluster 
% k: Rice k factor for LOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H_n : the channel (normalized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the steering vector function of a ULA of lambda/2 spacing along the y-axis (considering that theta = pi/2 
% in ISO Spherical coordinate system)
% the steering vector function of a ULA of lambda/2 spacing along the y-axis (considering that theta = pi/2 
% in ISO Spherical coordinate system)
a = @(phi,N) exp(-1j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);

% The NON-LOS contribution to the channel
[H_NLOS, Phi_AOD_NLOS, Phi_AOA_NLOS, Alpha_NLOS] = ...
    mm_wave_channel_v2_2D(Nt, Nr, Nc, Np, sig);
% The LOS contribution to the channel
phi_AOD_LOS = pi/2;
phi_AOA_LOS = pi/2;
At_LOS = a(phi_AOD_LOS, Nt); % Tx LOS steering vector
Ar_LOS = a(phi_AOA_LOS, Nr); % Rx LOS steering vector
H_LOS = sqrt(Nt*Nr)*Ar_LOS*At_LOS'; % LOS channel contribution (normalized)
H_n = sqrt(k/(k+1))*H_LOS + sqrt(1/(k+1))*H_NLOS; % channel

% This is just for verification prurposes
Phi_AOD = [Phi_AOD_NLOS(:); phi_AOD_LOS];
Phi_AOA = [Phi_AOA_NLOS(:); phi_AOA_LOS];
Alpha = [Alpha_NLOS(:); 1];             
end

