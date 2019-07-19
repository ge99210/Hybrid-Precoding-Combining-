%%%%%%%%%%%%%%%%%% Hybrid Precoding with HPSAGS %%%%%%%%%%%%%%%%
% Hybrid Precoding in mmWave channels with perfect 
% channel-state-information (PCSI) vs SNR
% Georgios K. Papageorgiou, 19/07/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
%close all;

N_iter = 10;

% Select between: full = 1 runs mdp or full = 0 does not run (faster)
full = 1;

tic;

Ns = 8; % number of transmitted streams
Nt = 64; % number of transmitter antennas
Mt = 8; % number of transmitter RF chains
Nr = 64; % number of receiver antennas
Mr = Mt; % number of receiver RF chains-not used

% SNR in dB (later converted to linear)
SNR_set = -30:5:0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPSAGS input parameters (fixed to)
Tmax = 100;
EPS = 1e-4;
M = 7;
mu_v(M) = 2.5;
for m = M:-1:2 
    mu_v(m-1) = mu_v(m)/2;
end
mu_v = flip(mu_v);
% the learning rate (most important parameter of SAGS) - 10 optimal for the precoding 
eta = 10; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel's characteristics
Nc = 5; 
Np = 3; 
sig = deg2rad(7.5); % this is the s.d. of the Laplace distribution for the scattering clusters

% Cost function for the steering vector of the ULA
a = @(phi,N) exp(-1j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);

% Initialization for the Matrix Decomposition Precoding
eps = 0.1;
delta = 0.1;
K_u = 100;

% Initialization
SE_Full_Dig = zeros(length(SNR_set),1); 
SE_Hybrid_omp = zeros(length(SNR_set),1);
SE_Hybrid_sags = zeros(length(SNR_set),1);
SE_Hybrid_md = zeros(length(SNR_set),1);

Time_omp = zeros(length(SNR_set),1);
Time_sags = zeros(length(SNR_set),1);
Time_md = zeros(length(SNR_set),1);

% Progress bar - comment while debugging
pbar=waitbar(0,'Please wait...','Name','Progress');

for SNR_index = 1:length(SNR_set)
    SNR = SNR_set(SNR_index); % in dB (SNR = 10*log10(rho))
    rho = 10^(SNR/10); % SNR (we assume that s_n = 1 for the noise
    
    Temp_se_full_dig = 0;
    Temp_se_sags = 0;
    Temp_se_omp = 0;
    Temp_se_md = 0;
        
    % Time initialization
    Tot_time_omp = 0;
    Tot_time_sags = 0;
    Tot_time_md = 0;
    
    for i=1:N_iter   
        
            % the mmWave channel
            [H, Phi_AOD, Phi_AOA, Alpha] =...
                        mm_wave_channel_v3_2D(Nt, Nr, Nc, Np, sig);
            % OMP's A matrices
            atv = Phi_AOD(:);
            arv = Phi_AOA(:);
            A_t = zeros(Nt,Nc*Np);
            A_r = zeros(Nr,Nc*Np);
            for m = 1:Nc*Np
               A_t(:,m) = a(atv(m), Nt);
               A_r(:,m) = a(arv(m), Nr);
            end   

        % Channel svd
        [U,S,V]=svd(H); 
        U1 = U(:,1:Ns);
        S1 = S(1:Ns,1:Ns);
        V1 = V(:, 1:Ns);

        % SE of a Fully Digital system (optimal)
        SE_full_dig = log2(det(eye(Ns)+(rho/Ns)*(S1.^2)));
        Temp_se_full_dig = Temp_se_full_dig + SE_full_dig;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % HP with OMP
        tic;
        [F_R_md, F_B_md] = hp_omp(V1, A_t, Ns, Mt);
        time_OMP = toc;
        R_omp = log2(det(eye(Nr)+(rho/Ns)*H*(F_R_md*F_B_md)*(F_B_md'*F_R_md')*H'));
        Temp_se_omp = Temp_se_omp + real(R_omp);
        Tot_time_omp = Tot_time_omp + time_OMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HP with HPSAGS
        Theta_0 = 2*pi*rand(Nt,Mt)-pi; % initialization
        tic;
        [F_R_sags, F_B_sags, cond_conv, k_prec]  =...
            hp_sags(V1, mu_v, eta, Theta_0, Tmax, EPS);
        time_sags = toc;
        F_sags = F_R_sags*F_B_sags;
        R_sags = log2(det(eye(Nr)+(rho/Ns)*H*(F_sags*F_sags')*H'));
        Temp_se_sags = Temp_se_sags + real(R_sags);
        Tot_time_sags = Tot_time_sags + time_sags;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % HP via Matrix Decomposition Precoding
        if full==1
        [U_F,S_F,V_F] = svd(V1);
        Ph_0 = angle(U_F*S_F);
        F_R_0 = exp(1j*Ph_0)/sqrt(Nt);
        tic;
        [F_R_md, F_B_md] = mdp(V1, F_R_0, eps, delta, K_u); % try the adaptive threshold
        time_md = toc;
        F_md = F_R_md*F_B_md;
        R_md = log2(det(eye(Nr)+(rho/Ns)*H*(F_md*F_md')*H'));
        Temp_se_md = Temp_se_md + R_md;
        Tot_time_md = Tot_time_md + time_md;   
        end
    end
       
    SE_Full_Dig(SNR_index) = real(Temp_se_full_dig)/N_iter;
    SE_Hybrid_omp(SNR_index) = real(Temp_se_omp)/N_iter;    
    SE_Hybrid_sags(SNR_index) = real(Temp_se_sags)/N_iter;   
    SE_Hybrid_md(SNR_index) = real(Temp_se_md)/N_iter;  
    
    Time_omp(SNR_index)= Tot_time_omp/N_iter;
    Time_sags(SNR_index)= Tot_time_sags/N_iter;
    Time_md(SNR_index)= Tot_time_md/N_iter;
    
    % Update waitbar and message
    fi=round(SNR_index*1000/length(SNR_set))/10;
    formatSpec = ' %1$3.1f %2$c';
    waitbar(fi/100,pbar,sprintf(formatSpec,fi,'%'));

end

close(pbar);
time_tot = toc/60; % in minutes

disp(['The maximum SE gain of SAGS over OMP is ', num2str(max(SE_Hybrid_sags-SE_Hybrid_omp)), ' bits/s/Hz']);
disp(['The maximum SE gain of SAGS over MD is ', num2str(max(SE_Hybrid_sags-SE_Hybrid_md)), ' bits/s/Hz']);
d_vec = SE_Hybrid_sags - SE_Hybrid_omp

f1 = figure(1);
movegui(f1,'west');
plot(SNR_set,SE_Full_Dig,'ko-', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,SE_Hybrid_omp,'rv--', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_md,'g^:', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_sags,'bo-.', 'Linewidth', 1.5,'MarkerSize',4);
hold off;
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
hl = legend('Optimal U-SVD','SSPOMP','MDP','HPSAGS','Location','Northwest');
set(hl, 'Fontsize', 12,'Interpreter','latex');
grid on;

f2 = figure(2);
movegui(f2,'east');
plot(SNR_set,Time_omp,'rv--', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,Time_md,'g^:', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,Time_sags,'bo-.', 'Linewidth', 1.5,'MarkerSize',4);
hold off;
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(gca, 'YScale', 'log');
set(get(gca,'YLabel'),'String','Aver. run time (sec)','Interpreter','latex');
h2 = legend('OMP', 'MD','SAGS','Location','Northwest');
set(h2, 'Fontsize', 12,'Interpreter','latex');
grid on;

%save('256x64MaMIMO_12RF_8Ns_Rayleigh_Imperfect_CSI_a=03_10MC_runs.mat','Capacity_Hybrid_EqPwrAlloc','Capacity_Hybrid_an','Capacity_Hybrid_ch','Capacity_Hybrid_Sohail');