function [W_R, W_B] = hc_omp(Y, W_mmse, A_r, M)
% Heath's OMP for hybrid precoding
% Georgios K. Papageorgiou 19/7/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ayach, O. El, Rajagopal, S., Abu-Surra, S., Pi, Z., & Heath, R. W. (2014). 
% Spatially sparse precoding in millimeter wave MIMO systems. IEEE Transactions on Wireless Communications, 13(3), 1499–1513. 
% https://doi.org/10.1109/TWC.2014.011714.130846
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_R = [];
W_res = W_mmse;
Y = eye(length(Y));
    for i=1:M
        Psi = A_r'*Y*W_res;
        D = diag(Psi*Psi');
        [~, k] = max(D);
        W_R = [W_R A_r(:,k)];
        W_B = (W_R'*Y*W_R)\W_R'*Y*W_mmse;
        W_res = (W_mmse-W_R*W_B)/norm(W_mmse-W_R*W_B,'fro');
    end
end

