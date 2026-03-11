function [H_AP_UE,HMean_AP_UE,HMean_RIS_UE,HMean_AP_RIS,H_AP_RIS,H_RIS_UE,Ngroup] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMeanWithoutPhase_AP_UE,HMeanWithoutPhase_AP_RIS,HMeanWithoutPhase_RIS_UE,groupRIS_size,N_H_RIS,N_V_RIS)
% Esta función genera canales Rician con correlación espacial entre AP-UE,
% RIS-AP y RIS-UE
% Argumentos de entrada:
% R_AP_UE                   = Matriz de correlación espacial entre el AP y el UE
%                             de dimensiones (L*N_AP,L*N_AP,L,K)
% R_AP_RIS1                 = Matriz de correlación espacial entre el AP y la RIS
%                             de dimensiones (N_AP,N_AP,L,S)
% R_AP_RIS2:                = Matriz de correlación espacial entre el AP y la RIS
%                             de dimensiones (N_RIS,N_RIS,L,S)
% R_RIS_UE                  = Matriz de correlación espacial entre la RIS y el UE
%                             de dimensiones (N_RIS,N_RIS,S,K)
% nbrOfRealizations         = Número de realizaciones del canal  
% L                         = Número de APs
% K                         = Número de usuarios
% S                         = Número de RIS
% N_AP                      = Número de antenas de cada AP
% N_RIS                     = Número de elementos de la RIS
% HMeanWithoutPhase_AP_UE   = Módulo de la componente de LoS del canal AP-UE
% HMeanWithoutPhase_AP_RIS  = Módulo de la componente de LoS del canal AP-RIS
% HMeanWithoutPhase_RIS_UE  = Módulo de la componente de LoS del canal RIS-UE
% groupRIS_size             = Tamaño de los grupos de una RIS
% N_H_RIS                   = Número de columnas de la RIS
% N_V_RIS                   = Número de filas de la RIS
% Argumentos de salida:
% H_AP_UE                   = Canal Rician generado entre el AP y el UE de
%                             dimensiones (L*N_AP,nbrOfRealizations,K)
% HMean_AP_UE               = Componente de LoS del canal AP-UE
%                             se usa agrupación se da el canal agrupado
% H_AP_RIS                  = Canal Rician generado entre el AP y la RIS
%                             se usa agrupación se da el canal agrupado
% H_RIS_UE                  = = Canal Rician generado entre el RIS y el UE
% Ngroup                    = Número de agrupaciones de una RIS
%% Generar realizaciones de canal

%----- AP-UE -----
% Generar canal Rician para AP-UE
M_AP_UE = L*N_AP;
H_AP_UE = zeros(M_AP_UE, nbrOfRealizations, K); % Canal resultante
W_AP_UE = (randn(M_AP_UE, nbrOfRealizations, K) + 1i * randn(M_AP_UE, nbrOfRealizations, K));   % Ruido gaussiano complejo

% Media del canal AP-UE
HMean_AP_UE=zeros(M_AP_UE,nbrOfRealizations,K); 
HMeanx_AP_UE=reshape(repmat(HMeanWithoutPhase_AP_UE,nbrOfRealizations,1),M_AP_UE,nbrOfRealizations,K);  % Se repite el canal tantas veces como realizaciones haya

% Fase aleatoria para componente LoS AP-UE
angles_AP_UE= -pi + 2*pi*rand(M_AP_UE,nbrOfRealizations,K);
phaseMatrix_AP_UE=exp(1i*angles_AP_UE);

% Canal Rician con correlación espacial AP-UE
for l = 1:L
    for k = 1:K
        HMean_AP_UE(:,:,k)= phaseMatrix_AP_UE(:,:,k).*HMeanx_AP_UE(:,:,k);  % Aplicar fase aleatoria
        Rsqrt = sqrtm(R_AP_UE(:,:,l,k));
        H_AP_UE((l-1)*N_AP+1:l*N_AP,:,k) = sqrt(0.5)*Rsqrt*W_AP_UE((l-1)*N_AP+1:l*N_AP,:,k) + HMean_AP_UE((l-1)*N_AP+1:l*N_AP,:,k);     % Se aplica un ruido con correlación al canal
       
    end
end

% ----- RIS-UE -----
% Generar canal Rician para RIS-UE
group_side = sqrt(groupRIS_size);
Ngroup = (N_H_RIS/group_side)*(N_V_RIS/group_side);

M_RIS_UE = S*N_RIS;
H_RIS_UE = zeros(M_RIS_UE, nbrOfRealizations, K); % Canal resultante
W_RIS_UE = (randn(M_RIS_UE, nbrOfRealizations, K) + 1i * randn(M_RIS_UE, nbrOfRealizations, K));   % Ruido gaussiano complejo

% Media del canal RIS-UE
HMean_RIS_UE=zeros(M_RIS_UE,nbrOfRealizations,K); 
HMeanx_RIS_UE=reshape(repmat(HMeanWithoutPhase_RIS_UE,nbrOfRealizations,1),M_RIS_UE,nbrOfRealizations,K);   % Se repite el canal tantas veces como realizaciones haya
 
% Fase aleatoria para componente LoS RIS-UE
angles_RIS_UE= -pi + 2*pi*rand(M_RIS_UE,nbrOfRealizations,K);
phaseMatrix_RIS_UE=exp(1i*angles_RIS_UE);

% Canal Rician con correlación espacial RIS-UE
for s = 1:S
    for k = 1:K
        
        HMean_RIS_UE(:,:,k)= phaseMatrix_RIS_UE(:,:,k).*HMeanx_RIS_UE(:,:,k);  % Aplicar fase aleatoria
        Rsqrt = sqrtm(R_RIS_UE(:,:,s,k));
        H_RIS_UE((s-1)*N_RIS+1:s*N_RIS,:,k) = sqrt(0.5)*Rsqrt*W_RIS_UE((s-1)*N_RIS+1:s*N_RIS,:,k) + HMean_RIS_UE((s-1)*N_RIS+1:s*N_RIS,:,k);    % Se aplica un ruido con correlación al canal        
    end
end

% ----- AP-RIS -----
% Generar canal Rician para AP-RIS
H_AP_RIS = zeros(L*N_AP, S*N_RIS, nbrOfRealizations);
W_AP_RIS = randn(L*N_AP, S*N_RIS, nbrOfRealizations) + 1i*randn(L*N_AP, S*N_RIS, nbrOfRealizations);

% Media del canal AP-RIS
HMean_AP_RIS = zeros(L*N_AP, S*N_RIS, nbrOfRealizations, 'single');
HMeanx_AP_RIS = reshape(repmat(HMeanWithoutPhase_AP_RIS, 1, nbrOfRealizations), L*N_AP, S*N_RIS, nbrOfRealizations);    % Se repite el canal tantas veces como realizaciones haya

% Fase aleatoria para componente LoS AP-RIS
angles_AP_RIS = -pi + 2*pi*rand(L*N_AP, S*N_RIS, nbrOfRealizations);
phaseMatrix_AP_RIS = exp(1i*angles_AP_RIS);

for l = 1:L
    for s = 1:S
        % Actualiza la media con fase aleatoria para el bloque (l,s)
        HMean_AP_RIS(:,(s-1)*N_RIS+1:s*N_RIS,:) = phaseMatrix_AP_RIS(:, (s-1)*N_RIS+1:s*N_RIS, :) .* HMeanx_AP_RIS(:, (s-1)*N_RIS+1:s*N_RIS, :);
        Rsqrt1 = sqrtm(R_AP_RIS1(:,:,l,s));
        Rsqrt2 = sqrtm(R_AP_RIS2(:,:,l,s));

        for t = 1:nbrOfRealizations
            % Multiplicación matricial sin squeeze, accediendo directo
            H_AP_RIS((l-1)*N_AP+1:l*N_AP, (s-1)*N_RIS+1:s*N_RIS, t) = sqrt(0.5)*Rsqrt1*W_AP_RIS((l-1)*N_AP+1:l*N_AP,(s-1)*N_RIS+1:s*N_RIS,t)*Rsqrt2 + HMean_AP_RIS((l-1)*N_AP+1:l*N_AP,(s-1)*N_RIS+1:s*N_RIS,t);    %Se añade al canal ruido con correlación espacial
        end
    end
end
