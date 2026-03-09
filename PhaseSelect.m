function [thetaMatrix] = PhaseSelect(nbrOfRealizations,H_AP_UE,Hhat_cascade,risAssignment,S,Ngroup,p,N_AP,L)

nbrOfIterations = 5;
thetaMatrix = exp(1i*(-pi + 2*pi*rand(Ngroup, S, nbrOfRealizations)));

for t = 1:nbrOfRealizations  % Por cada realización
    for s = 1:S             % Por cada RIS

        % Canal directo AP-UE, para realización t
        h_s = squeeze(H_AP_UE(:, t, risAssignment{s}));  
        for iter = 1:nbrOfIterations 
            for n = 1:Ngroup      % Por cada elemento n de la RIS s

                Hn = h_s + Hhat_cascade(:,(Ngroup * (s - 1) + 1:s*Ngroup),t,risAssignment{s})*thetaMatrix(:,s,t) - thetaMatrix(n,s,t)* Hhat_cascade(:,n + (Ngroup * (s - 1)),t,risAssignment{s})    ;%h_r(:,n)*h_t(n,:)
                bn = p*Hn;
                An = eye(N_AP*L) + p*(Hn*Hn') + p*Hhat_cascade(:,n + (Ngroup * (s - 1)),t,risAssignment{s})*(Hhat_cascade(:,n + (Ngroup * (s - 1)),t,risAssignment{s}))';

                thetaMatrix(n,s,t) = exp(-1i*angle(bn' *(An\Hhat_cascade(:,n + (Ngroup * (s - 1)),t,risAssignment{s}))));              
            end
        end
    end
end
