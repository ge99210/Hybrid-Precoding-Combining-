# Hybrid-Precoding-Combining-
Hybrid Processing (Precoding/Combining) for Rich and Poor Scattering Environments								

Cite this work as:

Georgios K. Papageorgiou, Mathini Sellathurai, Konstantinos Ntougias, Constantinos B. Papadias, "A Stochastic Optimization Approach 
to Hybrid Processing in Massive MIMO Systems," IEEE Wireless Communications Letters, 2019.

I. Precoding/Combining results in rich scattering environments (Rayleigh channels)

1. run hybrid_precoding_combining_Rayleigh_PCSI.m for comparison of HPSAGS with MDP under Perfect channel-state-information (PCSI) vs SNR

2. run hybrid_precoding_combining_Rayleigh_ICSI.m for comparison of HPSAGS with MDP under Imperfect channel-state-information (ICSI) vs SNR 

3. run SE_results_Rayleigh_vs_SNR.m to plot combined results (from 1&2)

4. run hybrid_precoding_combining_Rayleigh_vs_N_antennas.m for  for comparison of HPSAGS with MDP under Perfect channel-state-information (PCSI) vs the number of antennas (equal at the transmitter and receiver, i.e., N_t = N_r)

II. Precoding/Combining results in poor scattering environments - millimeter Wave (mmWave) channels

1. run hybrid_precoding_mmwave.m for comparison of HPSAGS with MDP and SSPOMP under PCSI vs SNR (precoding only)

2. run hybrid_precoding_combining_mmwave.m for comparison of HPSAGS with MDP and SSPOMP under PCSI vs SNR (precoding/combining)

3. run SE_results_mmWave_vs_SNR.m to plot combined results for Ns (# datastreams in 2)






