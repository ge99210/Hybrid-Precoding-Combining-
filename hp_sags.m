function  [F_R, F_B, conv_cond, k_t, Err2, Err]  = hp_sags(F_opt, mu_v, eta, Theta, Tmax, EPS)
% Hybrid Precoding with Stochastic Approximation with Gaussian Smoothing (SAGS)
% Georgios K. Papageorgiou 6/7/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT Parameters:
% F_opt: the optimal precoder of a fully digital system
% mu_v: the smoothing sequence in dicreasing order --> 0
% eta: the learning rate
% Theta: an initial guess of the phases
% Tmax: maximum number of iterations for SGDM (term. crit. 1)
% EPS:  accuracy of approximation(term. crit. 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT Parameters:
% F_R: Analogue precoder matrix
% F_B: Digital precoder matrix
% conv_cond: convergence condition for the Theta (should be dicreasing)
% k_t: SGDM steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization 
Ns = size(F_opt,2);
K = length(mu_v);
k = 0;
k_t = zeros(K,1);
conv_cond = zeros(1,K);
g = @(Theta) exp(1j*Theta)/sqrt(size(Theta,1));
Err=[];
    while k < K
        mu = mu_v(k+1);      
        t = 0; 
        cond2 = 1/EPS;
        F_R = g(Theta); 
        F_B = (F_R'*F_R)\F_R'*F_opt;
        
        while t < Tmax  && cond2>=EPS
          
          Theta_old = Theta;
          % Optimize the loads and weights  
          delta = normrnd(0,1,size(Theta));

          % Calculate the pertubated points
          Th_plus = Theta + mu*delta;
          Th_minus = Theta - mu*delta;
                               
          % Compute the exact gradient
          L_plus = gradient_g(Th_plus,F_opt,F_B,g);
          L_minus = gradient_g(Th_minus,F_opt,F_B,g);
                     
          % Compute the estimate of the gradient
          ksi = (1/2)*(L_plus+L_minus);  
          
          % Gradient descent with momentum for improved convergence
          Theta = Theta - eta*ksi; 
          
          t = t+1;
          %cond1 = norm(tau0*ksi,'fro'); 
          cond2 = norm(Theta-Theta_old,'fro')/norm(Theta_old,'fro');
          Err = [Err; cond2];
        end
       k = k+1;
       k_t(k) = t;
       conv_cond(k) = cond2;
       F_B = F_B*sqrt(Ns)/norm(F_R*F_B,'fro');
       F = g(Theta)*F_B;
       Err2(k) = norm(F-F_opt,'fro')/norm(F_opt,'fro');
    end
    
    F_R = g(Theta);
    F_B = (F_R'*F_R)\F_R'*F_opt;
    F_B = F_B*sqrt(Ns)/norm(F_R*F_B,'fro');
end

function G = gradient_g(X,V_opt,F_B,g)
          M = V_opt - g(X)*F_B;
          G = -2*real(1j*g(X).*(conj(M)*F_B.'));
end
