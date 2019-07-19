function [F_R, F_B] = hp_omp(F_opt, A_t, Ns, M)
% Heath's OMP for hybrid precoding
% Georgios K. Papageorgiou 19/7/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ayach, O. El, Rajagopal, S., Abu-Surra, S., Pi, Z., & Heath, R. W. (2014). 
% Spatially sparse precoding in millimeter wave MIMO systems. IEEE Transactions on Wireless Communications, 13(3), 1499–1513. 
% https://doi.org/10.1109/TWC.2014.011714.130846
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_R = [];
F_res = F_opt;
    for i=1:M
        Psi = A_t'*F_res;
        D = diag(Psi*Psi');
        [~, k] = max(D);
        F_R = [F_R A_t(:,k)];
        F_B = (F_R'*F_R)\F_R'*F_opt;
        F_res = (F_opt-F_R*F_B)/norm(F_opt-F_R*F_B,'fro');
    end
F_B = sqrt(Ns)*F_B/norm(F_R*F_B,'fro');
end

