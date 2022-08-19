function T_est = compute_T_values(psi, joint_psi, P, mu)

Nstates = length(psi(:,1));
m = length(psi(1,:));
T_est = zeros(Nstates,Nstates);
for ell = 2:m
    for s1 = 1:Nstates
        for s2 = 1:Nstates
            if P(s1,s2) == 0
                continue;
            end
            psi1 = joint_psi(s1,s2,ell);
            psi2 = psi(s1,ell-1);
           
            T_val = log2(psi1^(psi1/(mu(s1)*P(s1,s2)))/psi2^(psi2/mu(s1)));
            T_est(s1,s2) = T_est(s1,s2) + T_val/m;
        end
    end
end

end

