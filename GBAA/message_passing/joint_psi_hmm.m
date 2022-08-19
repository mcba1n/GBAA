function [joint_psi, log_joint_psi] = joint_psi_hmm(Nstates, m, log_post, P, P_Y_S, y, F, B) 

log_joint_psi = Inf*ones(Nstates,Nstates,m);
joint_psi = zeros(Nstates,Nstates,m);

log_P =  arrayfun(@eln, P);
log_P_Y_S =  arrayfun(@eln, P_Y_S);

for ell = 2:m
    for s1 = 1:Nstates
        for s2 = 1:Nstates
                log_gamma = log_P(s1, s2) + log_P_Y_S(y(ell),s2);
                log_joint_psi(s1,s2,ell) = F(s1,ell-1) + log_gamma + B(s2,ell);

                log_joint_psi(s1,s2,ell) = log_joint_psi(s1,s2,ell) - log_post;
                joint_psi(s1,s2,ell) = eexp(log_joint_psi(s1,s2,ell));
        end
    end
end


end

