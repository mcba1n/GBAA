function [B] = B_hmm(y, m, P, P_Y_S)

Nstates = size(P,1);
B = Inf*ones(Nstates,m);

log_P =  arrayfun(@eln, P);
log_P_Y_S = arrayfun(@eln,P_Y_S);

%% initialisation
for s = 1:Nstates
    B(s,m) = eln(1);
end

%% recursion
for ell = m-1:-1:1  
    for s = 1:Nstates   
        next_state = find(P(s,:) > 0);
        for i = 1:length(next_state)  
            log_gamma = log_P(s, next_state(i)) + log_P_Y_S(y(ell+1),next_state(i));
            B(s,ell) = elnsum(B(s,ell), elnproduct(log_gamma, B(next_state(i),ell+1)));
        end
    end
end

end

