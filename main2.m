% Vaciar espacio de trabajo y cerrar figuras
close all;
clear;

%% Setup de simulación
nbrOfSetups = 100;   % Número de escenarios
nbrOfRealizations = 1000;    % Número de realizaciones

L = 20;         % Número de APs
N_AP = 4;        % Antenas por AP
N_H_RIS = 8;    % Número de filas de la RIS
N_V_RIS = 8;    % Número de columnas de la RIS
N_RIS = N_V_RIS*N_H_RIS;     % Número de elementos de la RIS
K = 1;          % Número de UEs
tau_c = 20000;     % Longitud del bloque de coherencia
p = 100;         % Potencia de transmisión (mW)
fc = 3.5;         % Frecuencia (GHz)
LoS = 2;         % Linea de visión directa
% Desviación estándar angular en el modelo de dispersión local (en radianes)
ASD_varphi = deg2rad(15);  % angulo de azimut 
ASD_theta = deg2rad(15);  % angulo de elevación
groupRIS_size = 1;         % 1,4,16,64

% Arreglos 3D para guardar resultados por tipo de canal 
SE_PMMSE_DCC = zeros(K, nbrOfSetups, 6);  

%% Numero de RIS
%S_values = [0,5,10,20,50,75];
S_values = 1;
for s = 1:length(S_values)
    S = S_values(s);
    tau_p = K + S*(N_RIS/groupRIS_size+1);
    for n = 1:nbrOfSetups
        disp(['Setup ' num2str(n) '/' num2str(nbrOfSetups) ' asistido por ' num2str(S) ' RIS']);
    
        % Generar escenario
        [R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,pilotIndex,D,HMean_AP_UE, HMean_AP_RIS, HMean_RIS_UE, probLoS_AP_UE, probLoS_RIS_UE] = setup(L,K,N_AP,N_RIS,tau_p,n,ASD_varphi,ASD_theta,LoS,fc,S,N_H_RIS,N_V_RIS);
        
        % Asignacion de RIS
        if S == 0
            risAssignment = [];
        else
            risAssignment = assignRIS(probLoS_AP_UE, probLoS_RIS_UE);
        end
        
        % Generar canales
        [H_AP_UE,HMean_AP_UE,HMean_RIS_UE,HMean_AP_RIS,H_AP_RIS,H_RIS_UE,Ngroup] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,groupRIS_size,N_H_RIS,N_V_RIS);
        
        % Canal en cascada
        
        H_cascade = zeros(L*N_AP,S*N_RIS,nbrOfRealizations,K);
        for m = 1:S
            for k = 1:K
                h_ris_ue = H_RIS_UE((m-1)*N_RIS+1 : m*N_RIS, :, k); 
                for t = 1:nbrOfRealizations
                    for l = 1:L
                        h_ap_ris = H_AP_RIS((l-1)*N_AP+1 : l*N_AP, (m-1)*N_RIS+1 : m*N_RIS, t);
                        H_cascade((l-1)*N_AP+1 : l*N_AP,(m-1)*N_RIS+1 : m*N_RIS,t,k) = h_ap_ris .* (h_ris_ue(:, t).');
                    end
                end
            end
        end
        if(groupRIS_size > 1)
            H_cascade_grouped = RISGrouping(H_cascade, Ngroup, groupRIS_size, L,N_AP,S,nbrOfRealizations,K,N_RIS,N_V_RIS,N_H_RIS);
            thetaMatrix = zeros(N_RIS,S,nbrOfRealizations);
            [thetaMatrix_temp] = PhaseSelect(nbrOfRealizations,H_AP_UE,H_cascade_grouped,risAssignment,S,Ngroup,p,N_AP,L);
            g_V = sqrt(groupRIS_size);
            g_H = sqrt(groupRIS_size);
            
            for t = 1:nbrOfRealizations
                for m = 1:S
                    % 1. Reconstruir la grilla 2D del grupo (N_V_group x N_H_group)
                    theta_group_2D = reshape(thetaMatrix_temp(:,m,t), N_V_RIS/g_V, N_H_RIS/g_H);
                    
                    % 2. Expandir la grilla 2D al tamaño real de la RIS
                    theta_expanded_2D = repelem(theta_group_2D, g_V, g_H);
                    
                    % 3. Aplanar de vuelta a un vector 1D (columna por columna automáticamente)
                    thetaMatrix(:,m,t) = theta_expanded_2D(:); 
                end
            end
            % for t = 1:nbrOfRealizations
            %     for m = 1:S
            %         thetaMatrix_aux = repelem(reshape(thetaMatrix_temp(:,m,t),N_V_RIS/sqrt(groupRIS_size),N_H_RIS/sqrt(groupRIS_size)).',sqrt(groupRIS_size),sqrt(groupRIS_size));
            %         thetaMatrix(:,m,t) = reshape(thetaMatrix_aux.',N_RIS,1);
            %     end
            % end
        else
            [thetaMatrix] = PhaseSelect(nbrOfRealizations,H_AP_UE,H_cascade,risAssignment,S,N_RIS,p,N_AP,L);
        end
        
            % Canal agregado con las fases de las RISs configuradas
            HMean_eq_aux = zeros(L*N_AP,K, nbrOfRealizations);
            H_eq_aux = zeros(L*N_AP, K, nbrOfRealizations);
            h_reflected = zeros(L*N_AP,K,nbrOfRealizations);
            h_mean_reflected = zeros(L*N_AP,K,nbrOfRealizations);
            R_eq = zeros(N_AP, N_AP, L, K);
            for t = 1:nbrOfRealizations
                theta_aux = thetaMatrix(:,:,t);
                h_reflected(:,:,t) =  H_AP_RIS(:,:,t) * diag(theta_aux(:)) * squeeze(H_RIS_UE(:, t, :));  % dim: (1x L*N_AP)
                h_mean_reflected(:,:,t) = HMean_AP_RIS(:,:,t) * diag(theta_aux(:)) * squeeze(HMean_RIS_UE(:, t, :));
                H_eq_aux(:,:,t) = squeeze(H_AP_UE(:,t,:)) + h_reflected(:,:,t);
                HMean_eq_aux(:,:,t) = squeeze(HMean_AP_UE(:,t,:)) + h_mean_reflected(:,:,t);
                H_eq = permute(H_eq_aux,[1,3,2]);
                HMean_eq = permute(HMean_eq_aux,[1,3,2]);
                for k = 1:K
                    H_eq_aux1 = H_eq(:,t,k);
                    H_eq_aux2 = reshape(H_eq_aux1,N_AP,L);
                    HMean_eq_aux1 = HMean_eq(:,t,k);
                    HMean_eq_aux2 = reshape(HMean_eq_aux1,N_AP,L);
                    for l = 1:L
                        R_eq(:,:,l,k) = R_eq(:,:,l,k) + (H_eq_aux2(:,l)- HMean_eq_aux2(:,l)).*(H_eq_aux2(:,l)-HMean_eq_aux2(:,l))';
                    end
                end
            end
            R_eq = R_eq / nbrOfRealizations;
            % Estimar el canal agregado
            [Hhat_agregated,B_agregated,C_agregated] = channelEstimates(H_eq,HMean_eq,[],[],[],R_eq,[],nbrOfRealizations,L,K,N_AP,K,pilotIndex,p,risAssignment,S,[]);
            
            % % Calcular SE
            [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat_agregated,H_eq,D,B_agregated,C_agregated,tau_c,tau_p,nbrOfRealizations,N_AP,K,L,p,R_eq,pilotIndex);

        % 
        % 
        % Guardar resultados en la dimensión
        SE_PMMSE_DCC(:,n,s) = SE_P_MMSE;

        clear Hhat H_cascade Hhat_cascade Hhat_agregated H_eq B B_cascade B_agregated  C C_cascade C_agregated R_eq HMean_cascade Hhat_cascade_aux;
    end
end

% save('SE_noGroup_100Setups1RIS1UEnoPrelog4')
% save('SE_GroupOf4_100Setups1RIS1UEnoPrelog4')

%% Graficar resultados
% figure; hold on; box on;
% set(gca,'fontsize',16);
% 
% % P-MMSE 
% aux1 = SE_PMMSE_DCC(:,:,1); % 0 RIS
% aux2 = SE_PMMSE_DCC(:,:,2); % 5 RIS
% aux3 = SE_PMMSE_DCC(:,:,3); % 10 RIS
% aux4 = SE_PMMSE_DCC(:,:,4); % 20 RIS
% aux5 = SE_PMMSE_DCC(:,:,5); % 50 RIS
% aux6 = SE_PMMSE_DCC(:,:,6); % 75 RIS
% 
% 
% plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'k-', 'LineWidth', 2);
% plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'r-',  'LineWidth', 2);
% plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'g-', 'LineWidth', 2);
% plot(sort(aux4(:)), linspace(0,1,K*nbrOfSetups), 'b-', 'LineWidth', 2);
% plot(sort(aux5(:)), linspace(0,1,K*nbrOfSetups), 'm-',  'LineWidth', 2);
% plot(sort(aux6(:)), linspace(0,1,K*nbrOfSetups), 'y-', 'LineWidth', 2);
% 
% % % Ejes y leyenda
% xlabel('Spectral efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
% ylabel('CDF', 'Interpreter', 'Latex');
% %legend({'P-MMSE 1 RIS 1 UE'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
% legend({'P-MMSE 0 RIS', 'P-MMSE 5 RIS', 'P-MMSE 10 RIS', 'P-MMSE 20 RIS', 'P-MMSE 50 RIS','P-MMSE 75 RIS'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
% xlim([0 30]);
