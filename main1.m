% Vaciar espacio de trabajo y cerrar figuras
close all;
clear;

%% Setup de simulación
nbrOfSetups = 10;   % Número de escenarios
nbrOfRealizations = 100;    % Número de realizaciones

L = 20;         % Número de APs
N_AP = 4;        % Antenas por AP
N_H_RIS = 8;    % Número de filas de la RIS
N_V_RIS = 8;    % Número de columnas de la RIS
N_RIS = N_V_RIS*N_H_RIS;     % Número de elementos de la RIS
K = 2;          % Número de UEs
tau_c = 20000;     % Longitud del bloque de coherencia
p = 100;         % Potencia de transmisión (mW)
fc = 3.5;         % Frecuencia (GHz)
LoS = 2;         % Linea de visión directa
% Desviación estándar angular en el modelo de dispersión local (en radianes)
ASD_varphi = deg2rad(15);  % angulo de azimut 
ASD_theta = deg2rad(15);  % angulo de elevación
if (N_H_RIS == N_V_RIS && log2(N_H_RIS^2) == floor(log2(N_H_RIS^2)))
    groupRIS_size = 4;         % 1,4,16,64
else
    groupRIS_size_H = 1;
    groupRIS_size_V = 3;
    groupRIS_size = groupRIS_size_H*groupRIS_size_V;
end

% Arreglos 3D para guardar resultados por tipo de canal 
SE_PMMSE_DCC = zeros(K, nbrOfSetups, 6);  
%SE_MR_DIST   = zeros(K, nbrOfSetups, 6);

%% Numero de RIS
% S_values = [0,5,10,20,50,75];
S_values = [0,5,10];
for s = 1:length(S_values)
    S = S_values(s);
    tau_p = K + S*(N_RIS/groupRIS_size+1);
    for n = 1:nbrOfSetups
        %n = 4;
        disp(['Setup ' num2str(n) '/' num2str(nbrOfSetups) ' asistido por ' num2str(S) ' RIS']);
    
        % Generar escenario
        [R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,pilotIndex,D,HMean_AP_UE, HMean_AP_RIS, HMean_RIS_UE, probLoS_AP_UE, probLoS_RIS_UE] = setup(L,K,N_AP,N_RIS,tau_p,n,ASD_varphi,ASD_theta,LoS,fc,S,N_H_RIS,N_V_RIS);
        
        % Asignacion de RIS
        if S == 0
            risAssignment = [];
        else
            risAssignment = assignRIS(probLoS_AP_UE, probLoS_RIS_UE);
        end
        
        %Generar canales
        if (N_V_RIS == N_H_RIS   && log2(N_H_RIS^2) == floor(log2(N_H_RIS^2)))     % RIS cuadrada y de dimensión multiplo de 2
            [H_AP_UE,HMean_AP_UE,HMean_RIS_UE_real,HMean_RIS_UE,HMean_AP_RIS_real,HMean_AP_RIS,H_AP_RIS_grouped,H_AP_RIS,H_RIS_UE_grouped,H_RIS_UE,R_cascade,R_cascade_grouped_element,Ngroup,H_cascade] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,groupRIS_size, N_H_RIS,N_V_RIS);
        else                        % RIS rectangular
            [H_AP_UE,HMean_AP_UE,HMean_RIS_UE_real,HMean_RIS_UE,HMean_AP_RIS_real,HMean_AP_RIS,H_AP_RIS_grouped,H_AP_RIS,H_RIS_UE_grouped,H_RIS_UE,R_cascade,R_cascade_grouped_element,Ngroup,H_cascade] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,groupRIS_size, N_H_RIS,N_V_RIS,groupRIS_size_H,groupRIS_size_V);
        end
        H_cascade = zeros(L*N_AP,N_RIS*S,nbrOfRealizations,K);
        for k = 1:K
            for t = 1:nbrOfRealizations
                for l = 1:L
                    for a = (l-1)*N_AP+1:l*N_AP
                        for m = 1:S
                            for element = 1:N_RIS*S
                                H_cascade(a,element,t,k) = H_AP_RIS(a,element,t)*H_RIS_UE(element,t,k);
                            end
                        end
                    end
                end
            end
        end
        H_cascade = zeros(L*N_AP, N_RIS*S, nbrOfRealizations, K);
        HMean_cascade = zeros(L*N_AP, N_RIS*S, nbrOfRealizations, K);
        
        
        for k = 1:K
            for t = 1:nbrOfRealizations
                % Multiplicación elemento a elemento de las matrices completas
                % H_AP_RIS(:,:,t) es [L*N_AP x N_RIS*S]
                % H_RIS_UE(:,t,k).' es [1 x N_RIS*S] (transpuesto)
                % MATLAB expande automáticamente el vector fila para multiplicar cada antena
                H_cascade(:,:,t,k) = H_AP_RIS(:,:,t) .* (H_RIS_UE(:,t,k).');
                HMean_cascade(:,:,t,k) = HMean_AP_RIS_real(:,:,t) .* (HMean_RIS_UE_real(:,t,k).');
            end
        end
        if (groupRIS_size > 1)
            % Parámetros de la agrupación
            g_V = sqrt(groupRIS_size); % Elementos por fila en el grupo
            g_H = sqrt(groupRIS_size); % Elementos por columna en el grupo
            
            % Preparar la nueva matriz agrupada
            % Dimensiones: [Antenas_AP, Total_Grupos, Realizaciones, Usuarios]
            H_cascade_grouped = zeros(L*N_AP, Ngroup*S, nbrOfRealizations, K);
            HMean_cascade_grouped = zeros(L*N_AP, Ngroup*S, nbrOfRealizations, K);
            
            for k = 1:K
                for t = 1:nbrOfRealizations
                    for m = 1:S
                        % Extraer todos los elementos de ESTA RIS específica (m)
                        idx_ris = (m-1)*N_RIS + 1 : m*N_RIS;
                        
                        for a = 1:L*N_AP
                            % 1. Extraer el vector 1D de la cascada y pasarlo a grilla 2D (N_V_RIS x N_H_RIS)
                            % Nota: Usamos la misma lógica de reshape que tienes en channelGeneration
                            H_casc_2D = reshape(H_cascade(a, idx_ris, t, k), [N_V_RIS, N_H_RIS]).';
                            HMean_casc_2D = reshape(HMean_cascade(a, idx_ris, t, k), [N_V_RIS, N_H_RIS]).';
                            
                            % 2. Crear matriz temporal para el grupo
                            H_group_temp = zeros(N_V_RIS/g_V, N_H_RIS/g_H);
                            HMean_group_temp = zeros(N_V_RIS/g_V, N_H_RIS/g_H);
                            
                            % 3. Recorrer la grilla y sumar los bloques 2x2
                            for i = 1:(N_V_RIS/g_V)
                                for j = 1:(N_H_RIS/g_H)
                                    bloque = H_casc_2D((i-1)*g_V+1 : i*g_V, (j-1)*g_H+1 : j*g_H);
                                    bloque_mean = HMean_casc_2D((i-1)*g_V+1 : i*g_V, (j-1)*g_H+1 : j*g_H);
                                    H_group_temp(i, j) = sum(bloque, 'all'); % Suma los 4 elementos
                                    HMean_group_temp(i, j) = sum(bloque_mean, 'all');
                                end
                            end
                            
                            % 4. Volver a aplanar (vector 1D) y asignar a la matriz final
                            idx_group = (m-1)*Ngroup + 1 : m*Ngroup;
                            H_cascade_grouped(a, idx_group, t, k) = reshape(H_group_temp, [1, Ngroup]);
                            HMean_cascade_grouped(a, idx_group, t, k) = reshape(HMean_group_temp, [1, Ngroup]);
                        end
                    end
                end
            end
        end
        if (S>0)
            % Estimar canales
            if (groupRIS_size == 1)
                [Hhat,~,~,Hhat_cascade,~,~] = channelEstimates(H_AP_UE,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,H_cascade,R_AP_UE,R_cascade_grouped_element,nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,Ngroup,[]);
            else
                [Hhat,~,~,Hhat_cascade,~,~] = channelEstimates(H_AP_UE,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,H_cascade_grouped,R_AP_UE,R_cascade_grouped_element,nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,Ngroup,HMean_cascade_grouped);
            end
                % [Hhat,~,~,Hhat_cascade,~,~] = channelEstimates(H_AP_UE,HMean_AP_UE,HMean_AP_RIS_real,HMean_RIS_UE_real,H_AP_RIS,H_RIS_UE,R_AP_UE,R_cascade_grouped_element,nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,Ngroup,N_RIS);

            Hhat_cascade_aux = zeros(N_AP*L,N_RIS*S,nbrOfRealizations,K);
            for k = 1:K
                for l = 1:L
                    for m = 1:S
                        for t = 1:nbrOfRealizations
                            for a = (l-1)*N_AP+1:l*N_AP
                                if (N_V_RIS == N_H_RIS   && log2(N_H_RIS^2) == floor(log2(N_H_RIS^2)))
                                    H_aux_1 = repelem(reshape(Hhat_cascade(a,(m-1)*Ngroup+1:m*Ngroup,t,k),N_V_RIS/sqrt(groupRIS_size),N_H_RIS/sqrt(groupRIS_size)).',sqrt(groupRIS_size),sqrt(groupRIS_size));
                                    Hhat_cascade_aux(a,(m-1)*N_RIS+1:m*N_RIS,t,k) = reshape(H_aux_1.',N_RIS,1);
                                else
                                    H_aux_1 = repelem(reshape(Hhat_cascade(a,(m-1)*Ngroup+1:m*Ngroup,t,k),N_V_RIS/groupRIS_size_V,N_H_RIS/groupRIS_size_H).',groupRIS_size_V,groupRIS_size_H);
                                    Hhat_cascade_aux(a,(m-1)*N_RIS+1:m*N_RIS,t,k) = reshape(H_aux_1.',N_RIS,1);
                                end
                            end
                        end
                    end
                end
            end
            [thetaMatrix] = PhaseSelect(nbrOfRealizations,Hhat,Hhat_cascade_aux,risAssignment,S,N_RIS,p,N_AP,L);
           
            % Canal agregado con las fases de las RISs configuradas
            HMean_eq_aux = zeros(L*N_AP,K, nbrOfRealizations);
            H_eq_aux = zeros(L*N_AP, K, nbrOfRealizations);
            h_reflected = zeros(L*N_AP,K,nbrOfRealizations);
            h_mean_reflected = zeros(L*N_AP,K,nbrOfRealizations);
            R_eq = zeros(N_AP, N_AP, L, K);
        
            for t = 1:nbrOfRealizations
                theta_aux = thetaMatrix(:,:,t);
                h_reflected(:,:,t) =  H_AP_RIS(:,:,t) * diag(theta_aux(:)) * squeeze(H_RIS_UE(:, t, :));  % dim: (1x L*N_AP)
                h_mean_reflected(:,:,t) = HMean_AP_RIS_real(:,:,t) * diag(theta_aux(:)) * squeeze(HMean_RIS_UE_real(:, t, :));
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
            [Hhat_agregated,B_agregated,C_agregated] = channelEstimates(H_eq,HMean_eq,[],[],[],R_eq,[],nbrOfRealizations,L,K,N_AP,k,pilotIndex,p,risAssignment,S,[]);
            
            % % Calcular SE
            [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat_agregated,H_eq,D,B_agregated,C_agregated,tau_c,tau_p,nbrOfRealizations,N_AP,K,L,p,R_eq,pilotIndex);

        else
            %[Hhat,~,~,Hhat_cascade,~,~] = channelEstimates(H_AP_UE,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,H_cascade,R_AP_UE,R_cascade_grouped_element,nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,Ngroup,[]);
            [Hhat_agregated,B_agregated,C_agregated] = channelEstimates(H_AP_UE,HMean_AP_UE,[],[],[],R_AP_UE,[],nbrOfRealizations,L,K,N_AP,K,pilotIndex,p,risAssignment,S,[],[]);

             % % Calcular SE
            [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat_agregated,H_AP_UE,D,B_agregated,C_agregated,tau_c,tau_p,nbrOfRealizations,N_AP,K,L,p,R_AP_UE,pilotIndex);

        end
        % 
        % 
        % Guardar resultados en la dimensión
        SE_PMMSE_DCC(:,n,s) = SE_P_MMSE;
        %sum(SE_P_MMSE)
        %SE_MR_DIST(:,n,s)  = SE_MR_dist;

        clear Hhat H_cascade Hhat_cascade Hhat_agregated H_eq B B_cascade B_agregated  C C_cascade C_agregated R_eq;
    end
end
save('SE_GroupOf4_Theta15_Phi15_10RIS_10setups');
% save('SE_noGroup_100Setups1RIS1UEnoPrelogV7LOS0')
% save('SE_GroupOf4_1Setups1RIS1UEnoPrelogV7LOS0')
% save('SE_GroupOf4_100Setups1RIS1UEnoPrelogV7LOS0')
% save('SE_GroupOf64_100Setups1RIS1UEnoPrelogV7')

%% Graficar resultados
figure; hold on; box on;
set(gca,'fontsize',16);

% P-MMSE 
aux1 = SE_PMMSE_DCC(:,:,1); % 0 RIS
aux2 = SE_PMMSE_DCC(:,:,2); % 5 RIS
aux3 = SE_PMMSE_DCC(:,:,3); % 10 RIS
aux4 = SE_PMMSE_DCC(:,:,4); % 20 RIS
aux5 = SE_PMMSE_DCC(:,:,5); % 50 RIS
aux6 = SE_PMMSE_DCC(:,:,6); % 75 RIS


plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'k-', 'LineWidth', 2);
plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'r-',  'LineWidth', 2);
plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'g-', 'LineWidth', 2);
plot(sort(aux4(:)), linspace(0,1,K*nbrOfSetups), 'b-', 'LineWidth', 2);
plot(sort(aux5(:)), linspace(0,1,K*nbrOfSetups), 'm-',  'LineWidth', 2);
plot(sort(aux6(:)), linspace(0,1,K*nbrOfSetups), 'y-', 'LineWidth', 2);

% % Ejes y leyenda
xlabel('Spectral efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
ylabel('CDF', 'Interpreter', 'Latex');
%legend({'P-MMSE 1 RIS 1 UE'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
legend({'P-MMSE 0 RIS', 'P-MMSE 5 RIS', 'P-MMSE 10 RIS', 'P-MMSE 20 RIS', 'P-MMSE 50 RIS','P-MMSE 75 RIS'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
xlim([0 30]);
