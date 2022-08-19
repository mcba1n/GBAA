function [psi, log_psi] = psi_hmm(Nstates, m, log_post, F, B)

log_psi = Inf*ones(Nstates,m);
psi = zeros(Nstates,m);

for ell = 1:m
    for s = 1:Nstates
        log_psi(s,ell) = elnsum(log_psi(s,ell), elnproduct(F(s,ell), B(s, ell)));
    end

    log_psi(:,ell) = log_psi(:,ell) - log_post;
    for i = 1:length(log_psi(:,ell))
        psi(i,ell) = eexp(log_psi(i,ell));
    end
end

end

