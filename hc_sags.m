function  [W_R, W_B, conv_cond, k_t]  = hc_sags(W_opt, mu_v, eta, Theta, Tmax, EPS)
% Hybrid Combining with Stochastic Approximation with Gaussian Smoothing (SAGS)
% Georgios K. Papageorgiou 6/7/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT Parameters:
% W_opt: the optimal combiner of a fully digital system
% mu_v: the smoothing sequence in dicreasing order --> 0
% eta: the learning rate
% Theta: an initial guess of the phases
% Tmax: maximum number of iterations for SGDM(term. crit. 1)
% EPS:  accuracy of approximation(term. crit. 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT Parameters:
% W_R: Analogue combiner matrix
% W_B: Digital combiner matrix
% conv_cond: convergence condition for the gradient (should be dicreasing)
% k_t: SGDM steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization 
K = length(mu_v);
k = 0;
k_t = zeros(K,1);
conv_cond = zeros(1,K);
g = @(Theta) exp(1j*Theta)/sqrt(size(Theta,1));
    while k < K
        mu = mu_v(k+1);      
        t = 0; 
        cond2 = 1/EPS;
        W_R = g(Theta); 
        W_B = (W_R'*W_R)\W_R'*W_opt;
        
        while t < Tmax  && cond2>=EPS
          
          Theta_old = Theta;
          % Optimize the loads and weights  
          delta = normrnd(0,1,size(Theta));

          % Calculate the pertubated points
          Th_plus = Theta + mu*delta;
          Th_minus = Theta - mu*delta;
                               
          % Compute the exact gradient
          L_plus = gradient_g(Th_plus,W_opt,W_B,g);
          L_minus = gradient_g(Th_minus,W_opt,W_B,g);
                     
          % Compute the estimate of the gradient
          ksi = (1/2)*(L_plus+L_minus);  
          
          % Gradient descent with momentum for improved convergence
          Theta = Theta - eta*ksi; 
          
          t = t+1;
          %cond1 = norm(tau0*ksi,'fro'); 
          cond2 = norm(Theta-Theta_old,'fro')/norm(Theta_old,'fro');
        end
       k = k+1;
       k_t(k) = t;
       conv_cond(k) = cond2;
    end
    
    W_R = g(Theta);
    W_B = (W_R'*W_R)\W_R'*W_opt;
end

function G = gradient_g(X,V_opt,F_B,g)
          M = V_opt - g(X)*F_B;
          G = -2*real(1j*g(X).*(conj(M)*F_B.'));
end
