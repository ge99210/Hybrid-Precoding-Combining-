function [ F_R, F_B ] = alternate_approximation_hybrid_precoder( V1, F_R_0, eps , delta, K_u)
%  W. Ni, X. Dong, and W.-S. Lu, “Near-Optimal Hybrid Processing for Massive MIMO Systems via Matrix Decomposition,” 2015.
F_R = F_R_0;
F_B = (F_R_0'*F_R_0)\F_R_0'*V1;
Mt = size(F_B,1);
Nt = size(F_R,1);
Ns = size(F_B,2);
cond = norm(V1-F_R_0*F_B,'fro')/norm(V1,'fro');
options = optimoptions('quadprog','Display','off');%,'MaxIter',5000); % or lsqlin 
k = 0;
    while cond>eps && k<=K_u
        k = k+1;
        % Solve F_R 
        Q = V1-F_R*F_B;
        Ph = angle(F_R);
        D = zeros(Nt,Mt);
            for p=1:size(Q,1)
               q = (Q(p,:)).';
               G = ((1j/sqrt(Nt))*diag(exp(1j*Ph(p,:)))*F_B).';
               %[D(p,:),resnorm,residual,exitflag,output,lambda] = lsqlin(G,q,[],[],[],[],-delta*ones(Mt,1),delta*ones(Mt,1),[],options);
               H = real(G'*G); 
               f = -real(q'*G);
               [D(p,:), ~,~,~] = quadprog(H,f,[],[],[],[],-delta*ones(Mt,1),delta*ones(Mt,1),[],options);
            end
        F_R = F_R + (1j/sqrt(Nt))*D.*exp(1j*Ph);    
        F_B = (F_R'*F_R)\F_R'*V1;
        cond_old = cond;
        cond = norm(V1-F_R*F_B,'fro')/norm(V1,'fro');
        cd = abs(cond_old-eps);
        if cd>0.15
            delta = 1.25*delta;
        elseif cd<=0.15
            delta = 0.8*delta;
        end
        if delta<0.1
            delta = 0.1;
        elseif delta>0.5
            delta = 0.5;
        end
    end
F_B = F_B*sqrt(Ns)/norm(F_R*F_B,'fro');

end

