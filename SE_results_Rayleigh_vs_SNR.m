clear all;
clc;

load('256x64MIMO_10RF_8Ns_Rayleigh_PSCI_10MC_runs.mat');
SE_Hybrid_md_total_pcsi = SE_Hybrid_md_total;
SE_Hybrid_sags_total_pcsi = SE_Hybrid_sags_total;

load('256x64MIMO_10RF_8Ns_Rayleigh_ISCI_10MC_runs.mat');
SE_Full_Dig_pcsi = SE_Full_Dig_o;
SE_Full_Dig_icsi = SE_Full_Dig;
SE_Hybrid_md_total_icsi = SE_Hybrid_md_total;
SE_Hybrid_sags_total_icsi = SE_Hybrid_sags_total;

f1 = figure(1);
movegui(f1,'west');
plot(SNR_set,SE_Full_Dig_pcsi,'ko-', 'Linewidth', 1.5,'MarkerSize',5);
hold on;
plot(SNR_set,SE_Hybrid_md_total_pcsi,'ro:', 'Linewidth', 1.5,'MarkerSize',5);
plot(SNR_set,SE_Hybrid_sags_total_pcsi,'bo-.', 'Linewidth', 1.5,'MarkerSize',5);
plot(SNR_set,SE_Full_Dig_icsi,'ks-', 'Linewidth', 1.5,'MarkerSize',6);
plot(SNR_set,SE_Hybrid_md_total_icsi,'rs:', 'Linewidth', 1.5,'MarkerSize',6);
plot(SNR_set,SE_Hybrid_sags_total_icsi,'bs-.', 'Linewidth', 1.5,'MarkerSize',6);
hold off;
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
hl = legend('Optimal U-SVD PCSI','MDP PCSI','HPSAGS PCSI',...
                    'U-SVD ICSI','MDP ICSI','HPSAGS ICSI','Location','Northwest');
set(hl, 'Fontsize', 12,'Interpreter','latex');
grid on;
