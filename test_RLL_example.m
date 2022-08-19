clear all; close all; clc;
rng(1,'twister');
addpath('utils/');
addpath('utils/message_passing/');
addpath('utils/log_semiring/');

% RLL source state-space
A = [1, 1; 1, 0];
Nstates = size(A,1);

% IUD Markov source
[P_iud,mu_iud,H] = max_source_ent(A);

% Simulation parameters
p_vals = 0:0.05:0.5;
max_iters = 20;
eps = 0.001;
m = 1000000;
rates_vec = zeros(length(p_vals),max_iters);
optimised_rates_vec = zeros(1,length(p_vals));

% Run simulation
for p_idx = 1:length(p_vals)
    fprintf('Starting optimisation for p=%.2f:\n', p_vals(p_idx));
    P = P_iud;
    mu = mu_iud;
    for iter = 1:max_iters
        % Emission probabilities
        p = p_vals(p_idx);
        y_quant = [1,2];
        P_Y_S = [1-p, p;
                 p, 1-p];

        % Generate random observations
        [y,S] = hmmgenerate(m, P, P_Y_S');
        s_0 = 1;

        % Forward and backward probabilities
        [log_post,F_log] = F_hmm(y, m, P, P_Y_S, s_0);
        B_log = B_hmm(y, m, P, P_Y_S);

        % Smoothed probabilities
        [psi, log_psi] = psi_hmm(Nstates, m, log_post, F_log, B_log);
        [joint_psi,log_joint_psi] = joint_psi_hmm(Nstates, m, log_post, P, P_Y_S, y, F_log, B_log);

        % Compute information rate
        T_est = compute_T_values(psi, joint_psi, P, mu);
        I = compute_rate_from_T_values(T_est, P, mu);
        rates_vec(p_idx, iter) = I;
        fprintf('%.4f\n', I);

        % Update Markov source
        A_noisy = compute_noisy_adj(T_est, A);
        [P, mu, ~] = max_source_ent(A_noisy);

        % Check early stopping condition
        if iter > 1
            if abs(rates_vec(p_idx, iter)-rates_vec(p_idx, iter-1)) < eps
                break;
            end
        end
    end
    optimised_rates_vec(p_idx) = I;
    fprintf('Optimised rate: %.4f\n', I);
end

%% Theoretical bounds
H_max = log2(1+sqrt(5))-1;
C_b = @(p) 1 + p.*log2(p) + (1-p).*log2(1-p);
I_ub = C_b(p_vals);
I_ub(1) = 1;
I_lb = H_max*I_ub;

%% Plot the optimised rates with lower bound
plot(p_vals, rates_vec(:,1), '-k'), hold on, grid on;
plot(p_vals, I_lb, '--k'), hold on;
plot(p_vals, I_ub, '-.k');
ylim([0,1]);
xlim([0,0.5]);
legend('Optimised information rate', 'Theoretical lower bound', 'Theoretical upper bound')
xlabel('p');
ylabel('Information rate (bits/symbol)');
title('Information rates of BSC with RLL(0,1) source');

