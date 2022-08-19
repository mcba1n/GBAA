function [log_post,F] = F_hmm(y, m, P, P_Y_S, s_0)

Nstates = size(P,1);
F = Inf*ones(Nstates,m);

log_P =  arrayfun(@eln, P);
log_P_Y_S =  arrayfun(@eln, P_Y_S);

%% initialisation
for s = 1:Nstates
    log_gamma = log_P(s_0, s) + log_P_Y_S(y(1),s);
    F(s,1) = elnsum(F(s,1), log_gamma);
end

%% recursion
for ell = 2:m 
    for s = 1:Nstates
        prev_state = find(P(:,s) > 0);
        for i = 1:length(prev_state) 
            log_gamma = log_P(prev_state(i), s) + log_P_Y_S(y(ell), s);
            F(s,ell) = elnsum(F(s,ell), elnproduct(log_gamma, F(prev_state(i),ell-1)));  
        end
    end
end

%% termination
log_post = Inf;
for s = 1:Nstates
    log_post = elnsum(log_post, F(s,m));
end

end

