function [H_AP_UE,HMean_AP_UE,HMean_RIS_UE,HMean_RIS_UE_grouped,HMean_AP_RIS,HMean_AP_RIS_grouped,H_AP_RIS_grouped,H_AP_RIS,H_RIS_UE_grouped,H_RIS_UE,R_cascade,R_cascade_grouped_element,Ngroup,H_cascade] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMeanWithoutPhase_AP_UE,HMeanWithoutPhase_AP_RIS,HMeanWithoutPhase_RIS_UE,groupRIS_size,N_H_RIS,N_V_RIS,groupRIS_size_H,groupRIS_size_V)
% Esta función genera canales Rician con correlación espacial entre AP-UE,
% RIS-AP y RIS-UE y en caso de agrupar los elementos de la RIS genera los
% canales equivalentes de la agrupación
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
% HMean_RIS_UE_grouped      = Componente de LoS del canal RIS-UE
% HMean_AP_RIS_grouped      = Componente de LoS del canal AP-RIS
% H_AP_RIS_grouped          = Canal Rician generado entre el AP y la RIS si
%                             se usa agrupación se da el canal agrupado
% H_AP_RIS                  = Canal Rician generado entre el AP y la RIS
% H_RIS_UE_grouped          = Canal Rician generado entre la RIS y el UE si
%                             se usa agrupación se da el canal agrupado
% H_RIS_UE                  = = Canal Rician generado entre el RIS y el UE
% R_cascade                 = Correlación del canal en cascada
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
if (N_H_RIS == N_V_RIS && log2(N_H_RIS*N_V_RIS) == floor(log2(N_H_RIS*N_V_RIS)))
    group_side = sqrt(groupRIS_size);
    if (group_side > 1)
        Ngroup = (N_H_RIS / group_side)^2;
    else
        Ngroup = N_H_RIS^2;
    end
elseif(groupRIS_size_H>0 && groupRIS_size_V>0)
    Ngroup = (N_H_RIS/groupRIS_size_H)*(N_V_RIS/groupRIS_size_V);
else
    Ngroup = N_H_RIS*N_V_RIS;
end

M_RIS_UE = S*N_RIS;
H_RIS_UE = zeros(M_RIS_UE, nbrOfRealizations, K); % Canal resultante
H_RIS_UE_grouped = zeros(S*Ngroup, nbrOfRealizations, K);
W_RIS_UE = (randn(M_RIS_UE, nbrOfRealizations, K) + 1i * randn(M_RIS_UE, nbrOfRealizations, K));   % Ruido gaussiano complejo

% Media del canal RIS-UE
HMean_RIS_UE=zeros(M_RIS_UE,nbrOfRealizations,K); 
HMean_RIS_UE_grouped=zeros(S*Ngroup,nbrOfRealizations,K); 
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
H_AP_RIS_grouped = zeros(L*N_AP, S*Ngroup, nbrOfRealizations);
W_AP_RIS = randn(L*N_AP, S*N_RIS, nbrOfRealizations) + 1i*randn(L*N_AP, S*N_RIS, nbrOfRealizations);

% Media del canal AP-RIS
HMean_AP_RIS = zeros(L*N_AP, S*N_RIS, nbrOfRealizations, 'single');
HMean_AP_RIS_grouped = zeros(L*N_AP, S*Ngroup, nbrOfRealizations, 'single');
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
% %% Calculo de R_cascade
% 
% % Inicializamos 
% H_cascade_aux = zeros(M_AP_UE, K, nbrOfRealizations);
% R_cascade = zeros(N_AP, N_AP, L, K);
% R_cascade_grouped_element = zeros(N_AP, N_AP, L, K, S*Ngroup);
% h_reflected = zeros(L*N_AP,K,nbrOfRealizations);
% % 
% for t = 1:nbrOfRealizations  % Por cada realización
%     h_reflected(:,:,t) =  H_AP_RIS(:,:,t) * squeeze(H_RIS_UE(:, t, :));  % dim: (1x L*N_AP)
%     H_cascade_aux(:,:,t) = squeeze(H_AP_UE(:,t,:)) + h_reflected(:,:,t);
%     H_cascade = permute(H_cascade_aux,[1,3,2]);
%     for k = 1:K
%         H_eq_aux1 = H_cascade(:,t,k);
%         H_eq_aux2 = reshape(H_eq_aux1,N_AP,L);
%         for l = 1:L
%             R_cascade(:,:,l,k) = R_cascade(:,:,l,k) + (H_eq_aux2(:,l)- HMean_AP_UE((l-1)*N_AP+1 : l*N_AP, t, k)).*(H_eq_aux2(:,l)-HMean_AP_UE((l-1)*N_AP+1 : l*N_AP, t, k))';
% 
%             % Covarianza INDIVIDUAL de cada grupo de la RIS
%             for s = 1:S
%                 for j = 1:Ngroup
% 
%                     % Extraer el canal de ese grupo específico y su media
%                     h_g = H_AP_RIS_grouped((l-1)*N_AP+1:l*N_AP, j + Ngroup * (s - 1), t) * H_RIS_UE_grouped(j + Ngroup * (s - 1), t, k);
%                     h_g_mean = HMean_AP_RIS_grouped((l-1)*N_AP+1:l*N_AP, j + Ngroup * (s - 1), t) * HMean_RIS_UE_grouped(j + Ngroup * (s - 1), t, k);
% 
%                     % Calcular la desviación y sumar al acumulador de covarianza
%                     R_cascade_grouped_element(:,:,l,k,j + Ngroup * (s - 1)) = R_cascade_grouped_element(:,:,l,k,j + Ngroup * (s - 1)) + ((h_g-h_g_mean)*(h_g-h_g_mean)');
%                 end
%             end
%         end
%     end
% end
% R_cascade = R_cascade / nbrOfRealizations;
% R_cascade_grouped_element = R_cascade_grouped_element / nbrOfRealizations;

%% Calculo de R_cascade físico (Suma de productos)

% 1. Crear mapa de pertenencia a grupos (1D) para la RIS
group_map = zeros(N_V_RIS, N_H_RIS);
idx = 1;
if (N_H_RIS == N_V_RIS && log2(N_H_RIS*N_V_RIS) == floor(log2(N_H_RIS*N_V_RIS)))
    g_V = sqrt(groupRIS_size);
    g_H = sqrt(groupRIS_size);
else
    g_V = groupRIS_size_V;
    g_H = groupRIS_size_H;
end

for i = 1 : N_V_RIS/g_V
    for j = 1 : N_H_RIS/g_H
        group_map((i-1)*g_V + 1 : i*g_V, (j-1)*g_H + 1 : j*g_H) = idx;
        idx = idx + 1;
    end
end
group_map_1d_aux = group_map.'; % Vector que indica a qué grupo pertenece cada elemento (1 a N_RIS)
group_map_1d = group_map_1d_aux(:);
% 2. Inicializamos
H_cascade_aux = zeros(M_AP_UE, K, nbrOfRealizations);
R_cascade = zeros(N_AP, N_AP, L, K);
R_cascade_grouped_element = zeros(N_AP, N_AP, L, K, S*Ngroup);
h_reflected = zeros(L*N_AP,K,nbrOfRealizations);

% 3. Calcular covarianzas
for t = 1:nbrOfRealizations  
    h_reflected(:,:,t) =  H_AP_RIS(:,:,t) * squeeze(H_RIS_UE(:, t, :));  
    H_cascade_aux(:,:,t) = squeeze(H_AP_UE(:,t,:)) + h_reflected(:,:,t);
    H_cascade = permute(H_cascade_aux,[1,3,2]);
    for k = 1:K
        H_eq_aux1 = H_cascade(:,t,k);
        H_eq_aux2 = reshape(H_eq_aux1,N_AP,L);
        for l = 1:L
            % Covarianza del canal total
            R_cascade(:,:,l,k) = R_cascade(:,:,l,k) + (H_eq_aux2(:,l)- HMean_AP_UE((l-1)*N_AP+1 : l*N_AP, t, k)).*(H_eq_aux2(:,l)-HMean_AP_UE((l-1)*N_AP+1 : l*N_AP, t, k))';

            % Covarianza FISICA de cada grupo de la RIS
            for s = 1:S
                for j = 1:Ngroup
                    % Obtener los índices exactos de los elementos del grupo j
                    idx_elems = find(group_map_1d == j) + (s-1)*N_RIS;

                    % h_g real es la SUMA de los productos (AP-RIS * RIS-UE)
                    h_g = sum( H_AP_RIS((l-1)*N_AP+1:l*N_AP, idx_elems, t) .* (H_RIS_UE(idx_elems, t, k).'), 2);
                    h_g_mean = sum( HMean_AP_RIS((l-1)*N_AP+1:l*N_AP, idx_elems, t) .* (HMean_RIS_UE(idx_elems, t, k).'), 2);

                    % Calcular la desviación y sumar al acumulador
                    R_cascade_grouped_element(:,:,l,k,j + Ngroup * (s - 1)) = R_cascade_grouped_element(:,:,l,k,j + Ngroup * (s - 1)) + ((h_g-h_g_mean)*(h_g-h_g_mean)');
                end
            end
        end
    end
end
R_cascade = R_cascade / nbrOfRealizations;
R_cascade_grouped_element = R_cascade_grouped_element / nbrOfRealizations;
