% test mamimo optimization 
%clc;
clear all;
%close all;

N_iter = 10;
% rng('default');

full = 0;

tic;

Ns = 8; % number of transmitted streams - 7 
Nt = 256; % number of transmitter antennas - 64
Mt = 10; % number of transmitter RF chains - 10
Nr = 64; % number of receiver antennas - 64 gain 3.9 b/s/Hz
Mr = Mt; % number of receiver RF chains-not used

% SNR in dB (later converted to linear)
SNR_set = -30:5:0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAGS input parameters (fixed to)
Tmax = 100;
EPS = 1e-4;
K = 7;
mu_v(K) = 2.5; 
for n = K:-1:2 
    mu_v(n-1) = mu_v(n)/2;
end
mu_v = flip(mu_v);
% the learning rate (most important parameter of SAGS) - 10 optimal for the
% precoding 600 for the combining
eta_p = 10;
eta_c = 600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function for the steering vector of the ULA
a = @(phi,N) exp(-1j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);

% Initialization for the Matrix Decomposition 
eps = 0.1;
delta = 0.1;
K_u = 100;

% Initialization
SE_Full_Dig = zeros(length(SNR_set),1); 
SE_Hybrid_omp_total = zeros(length(SNR_set),1);
SE_Hybrid_sags_total = zeros(length(SNR_set),1);
SE_Hybrid_md_total = zeros(length(SNR_set),1);

Time_omp = zeros(length(SNR_set),1);
Time_sags = zeros(length(SNR_set),1);
Time_md = zeros(length(SNR_set),1);

% Progress bar - comment while debugging
pbar=waitbar(0,'Please wait...','Name','Progress');

for SNR_index = 1:length(SNR_set)
    SNR = SNR_set(SNR_index); % in dB (SNR = 10*log10(rho))
    rho = 10^(SNR/10); % SNR (we assume that s_n = 1 for the noise)
    
    Temp_se_full_dig = 0;
    Temp_total_sags = 0;
    Temp_total_md = 0;
    Temp_se_ben = 0; % remove when done
    
    % Time initialization
    Tot_time_sags = 0;
    Tot_time_md = 0;

    for i=1:N_iter   
        
        H = (1/sqrt(2))*(randn(Nr,Nt)+1j*randn(Nr,Nt));
        
        % Channel svd
        [U,S,V]=svd(H); 
        U1 = U(:,1:Ns);
        S1 = S(1:Ns,1:Ns);
        V1 = V(:, 1:Ns);
        
%%%%%%%%%%%%%%%%%%%%%%%%% FULLY DIGITAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SE of a Fully Digital system (optimal)
        SE_full_dig = log2(det(eye(Ns)+(rho/Ns)*(U1'*H*V1)*(U1'*H*V1)'));
        Temp_se_full_dig = Temp_se_full_dig + SE_full_dig;      
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Theta_0_t = 2*pi*rand(Nt,Mt)-pi; % initialization
        % Solve the hybrid precoding with SAGS
        tic;
        [F_R_sags, F_B_sags, cond_prec, k_prec,Theta]  =...
            hp_sags(V1, mu_v, eta_p, Theta_0_t, Tmax, EPS);
        time_SAGSp = toc;
        F_sags = F_R_sags*F_B_sags;
        %F_sags = V1; % optimize the combiner's parameters 
        
        % The matrices for the combining       
        A_mmse_sags = (rho/Ns)*H*(F_sags*F_sags')*H'+eye(Nr);
        W_mmse_sags = A_mmse_sags\((sqrt(rho)/Ns)*H*F_sags);
      
        Theta_0_r = 2*pi*rand(Nr,Mr)-pi;  % initialization
        tic;
        [W_R_sags, W_B_sags, conv_comb, k_comb]  = ...
            hc_sags(W_mmse_sags, mu_v, eta_c, Theta_0_r, Tmax, EPS);
        time_SAGSc = toc;
        Tot_time_sags = Tot_time_sags + (time_SAGSp + time_SAGSc); 
        
        R_sags = inv(W_B_sags'*(W_R_sags'*W_R_sags)*W_B_sags);
        SE_total_sags = ...
            log2(det(eye(Ns)+(rho/Ns)*R_sags*(W_B_sags'*W_R_sags'*(H*F_sags)*F_sags')*(W_B_sags'*W_R_sags'*(H*F_sags)*F_sags')'));
        Temp_total_sags = Temp_total_sags + real(SE_total_sags);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
        
        if full==1
        % Matrix Decomposition 
        % precoding
        [U_F,S_F,V_F] = svd(V1);
        Ph_0_t = angle(U_F*S_F);
        F_R_0 = exp(1j*Ph_0_t)/sqrt(Nt);
        tic;
        [F_R_md, F_B_md] = ...
            alternate_approximation_hybrid_precoder(V1, F_R_0, eps, delta, K_u); 
        time_md_p = toc;
        
        % combining
        F_md = F_R_md*F_B_md;
        % The matrices for the combining       
        A_mmse_md = (rho/Ns)*H*(F_md*F_md')*H'+eye(Nr);
        W_mmse_md = A_mmse_md\((sqrt(rho)/Ns)*H*F_md);
        [U_W,S_W,V_W] = svd(W_mmse_md);
        Ph_0_r = angle(U_W*S_W);
        W_R_0 = exp(1j*Ph_0_r)/sqrt(Nr);
        tic;
        [W_R_md, W_B_md] = ...
            alternate_approximation_hybrid_precoder(W_mmse_md, W_R_0, eps, delta, K_u); 
        time_md_c = toc;
        Tot_time_md = Tot_time_md + (time_md_p+time_md_c);
        
        R_md = inv(W_B_md'*(W_R_md'*W_R_md)*W_B_md);
        
        SE_total_md = ...
            log2(det(eye(Ns)+(rho/Ns)*R_md*((W_B_md'*W_R_md')*H*(F_md*F_md'))*((W_B_md'*W_R_md')*H*(F_md*F_md'))'));
        Temp_total_md = Temp_total_md + real(SE_total_md);
        end
    end
       
    SE_Full_Dig(SNR_index) = real(Temp_se_full_dig)/N_iter;  
    SE_Hybrid_sags_total(SNR_index) = real(Temp_total_sags)/N_iter;   
    SE_Hybrid_md_total(SNR_index) = real(Temp_total_md)/N_iter;  
    
    Time_sags(SNR_index)= Tot_time_sags/N_iter;
    Time_md(SNR_index)= Tot_time_md/N_iter;
    
    % Update waitbar and message
    fi=round(SNR_index*1000/length(SNR_set))/10;
    formatSpec = ' %1$3.1f %2$c';
    waitbar(fi/100,pbar,sprintf(formatSpec,fi,'%'));

end

close(pbar);
time_tot = toc/60; % in minutes

disp(['The maximum SE gain of MD over SAGS is ', num2str(max(SE_Hybrid_md_total-SE_Hybrid_sags_total)), ' bits/s/Hz']);
d_vec2 = SE_Full_Dig - SE_Hybrid_sags_total 

f1 = figure(1);
movegui(f1,'west');
plot(SNR_set,SE_Full_Dig,'ko-', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,SE_Hybrid_md_total,'r^:', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_sags_total,'bo--', 'Linewidth', 1.5,'MarkerSize',4);
hold off;
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
hl = legend('Fully Digital','MD','SAGS','Location','Northwest');
set(hl, 'Fontsize', 12,'Interpreter','latex');
grid on;

%save('256x64MIMO_10RF_8Ns_Rayleigh_PSCI_10MC_runs.mat','SNR_set','SE_Full_Dig','SE_Hybrid_md_total', 'SE_Hybrid_sags_total');