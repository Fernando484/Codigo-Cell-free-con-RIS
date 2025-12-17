function [thetaMatrix] = PhaseSelect(nbrOfRealizations,thetaMatrix,H_AP_UE,Hhat_cascade,risAssignment,S,Ngroup,p,N_AP,L)

for t = 1:nbrOfRealizations  % Por cada realización
    for s = 1:S             % Por cada RIS
        
        % Canal directo AP-UE, para realización t
        h_s = squeeze(H_AP_UE(:, t, risAssignment{s}));  % dim: (M_AP_UE x 1)

        for n = 1:Ngroup      % Por cada elemento n de la RIS s

            Hn = h_s + Hhat_cascade(:,(Ngroup * (s - 1) + 1),t,risAssignment{s}) - thetaMatrix(n,s)* Hhat_cascade(:,n + (Ngroup * (s - 1)),t,risAssignment{s})    ;%h_r(:,n)*h_t(n,:)
            bn = p*Hn'*Hhat_cascade(:,(Ngroup * (s - 1) + 1),t,risAssignment{s});
            An = eye(N_AP*L) + p*(Hn*Hn') + p*Hhat_cascade(:,(Ngroup * (s - 1) + 1),t,risAssignment{s})*(Hhat_cascade(:,(Ngroup * (s - 1) + 1),t,risAssignment{s}))';

            thetaMatrix(n,s) = exp(-1i*angle(bn' *(An(:,n)\Hhat_cascade(:,(Ngroup * (s - 1) + 1),t,risAssignment{s}))));              
        end        
    end
end