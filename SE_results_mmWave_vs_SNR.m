clear all;
clc;

load('64x64MIMO_10RF_10Ns_mmWave_30MC_runs.mat');
SE_Full_Dig_10 = SE_Full_Dig;
SE_Hybrid_omp_total_10 = SE_Hybrid_omp_total;
SE_Hybrid_md_total_10 = SE_Hybrid_md_total;
SE_Hybrid_sags_total_10 = SE_Hybrid_sags_total;

load('64x64MIMO_10RF_5Ns_mmWave_30MC_runs.mat');
SE_Full_Dig_7 = SE_Full_Dig;
SE_Hybrid_omp_total_7 = SE_Hybrid_omp_total;
SE_Hybrid_md_total_7 = SE_Hybrid_md_total;
SE_Hybrid_sags_total_7 = SE_Hybrid_sags_total;

f1 = figure(1);
subplot(2,1,1);
plot(SNR_set,SE_Full_Dig_7,'ko-', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,SE_Hybrid_omp_total_7,'gv--', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_md_total_7,'r^:', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_sags_total_7,'bo-.', 'Linewidth', 1.5,'MarkerSize',4);
%set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
h2 = legend('Optimal U-SVD','SSPOMP','MDP','HPSAGS','Location','Northwest');
set(h2, 'Fontsize', 12,'Interpreter','latex');
grid on;
subplot(2,1,2);
plot(SNR_set,SE_Full_Dig_10,'ko-', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,SE_Hybrid_omp_total_10,'gv--', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_md_total_10,'r^:', 'Linewidth', 1.5,'MarkerSize',4);
plot(SNR_set,SE_Hybrid_sags_total_10,'bo-.', 'Linewidth', 1.5,'MarkerSize',4);
hold off;
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
% hl = legend('Optimal U-SVD','SSPOMP','MDP','HPSAGS','Location','Northwest');
% set(hl, 'Fontsize', 12,'Interpreter','latex');
grid on;
