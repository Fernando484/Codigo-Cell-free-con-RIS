% Vaciar espacio de trabajo y cerrar figuras
close all;
clear;

%% Setup de simulación
nbrOfSetups = 20;   % Número de escenarios
nbrOfRealizations = 100;    % Número de realizaciones

L = 20;         % Número de APs
N_AP = 4;        % Antenas por AP
N_H_RIS = 8;    % Número de filas de la RIS
N_V_RIS = N_H_RIS;    % Número de columnas de la RIS
N_RIS = N_V_RIS*N_H_RIS;     % Número de elementos de la RIS
K = 10;          % Número de UEs
tau_c = 20000;     % Longitud del bloque de coherencia
%tau_p = 10;      % Longitud del piloto
p = 100;         % Potencia de transmisión (mW)
fc = 3.5;         % Frecuencia (GHz)
LoS = 2;         % Linea de visión directa
% Desviación estándar angular en el modelo de dispersión local (en radianes)
ASD_varphi = deg2rad(15);  % angulo de azimut 
%ASD_theta = deg2rad(15);  % angulo de elevación
groupRIS_size = 1;         % 1,4,16

% Arreglos 3D para guardar resultados por tipo de canal 
SE_PMMSE_DCC = zeros(K, nbrOfSetups, 6);  
%SE_MR_DIST   = zeros(K, nbrOfSetups, 6);

%% Numero de RIS
S_values = [0,5,10,20,50,75];
%S_values = 5;
for s = 1:length(S_values)
    S = S_values(s);
    tau_p = K;
    for n = 1:nbrOfSetups
        disp(['Setup ' num2str(n) '/' num2str(nbrOfSetups) ' asistido por ' num2str(S) ' RIS']);
    
        % Generar escenario
        [R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,pilotIndex,D,HMean_AP_UE, HMean_AP_RIS, HMean_RIS_UE, probLoS_AP_UE, probLoS_RIS_UE] = setup(L,K,N_AP,N_RIS,tau_p,n,ASD_varphi,LoS,fc,S,N_H_RIS,N_V_RIS);
        
        % Asignacion de RIS
        if S == 0
            risAssignment = [];
        else
            risAssignment = assignRIS(probLoS_AP_UE, probLoS_RIS_UE);
        end
        
        %Generar canales
        [H_AP_UE,HMean_AP_UE,HMean_RIS_UE,HMean_AP_RIS,H_AP_RIS_grouped,H_AP_RIS,H_RIS_UE_grouped,H_RIS_UE,R_cascade,Ngroup,H_cascade] = channelGeneration(R_AP_UE,R_AP_RIS1,R_AP_RIS2,R_RIS_UE,nbrOfRealizations,L,K,S,N_AP,N_RIS,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,groupRIS_size, N_H_RIS);
        if (S>0)
            thetaMatrix = zeros(Ngroup,S);
            % Fase aleatoria
            for t = 1:nbrOfRealizations
                for l = 1:S
                    thetaMatrix(:,l) = exp(1i*2*pi*rand(Ngroup,1));
                end
            end
            thetaMatrix = repelem(thetaMatrix,groupRIS_size,1);
        
            % Canal agregado con las fases de las RISs configaudas
            H_eq_aux = zeros(L*N_AP, K, nbrOfRealizations);
            R_eq = zeros(N_AP, N_AP, L, K);
            for t = 1:nbrOfRealizations  % Por cada realización
                h_reflected(:,:,t) =  H_AP_RIS(:,:,t) * diag(thetaMatrix(:)) * squeeze(H_RIS_UE(:, t, :));  % dim: (1x L*N_AP)
                H_eq_aux(:,:,t) = squeeze(H_AP_UE(:,t,:)) + h_reflected(:,:,t);
                H_eq = permute(H_eq_aux,[1,3,2]);
                for k = 1:K
                    H_eq_aux1 = H_eq(:,t,k);
                    H_eq_aux2 = reshape(H_eq_aux1,N_AP,L);
                    for l = 1:L
                        R_eq(:,:,l,k) = R_eq(:,:,l,k) + (H_eq_aux2(:,l)- HMean_AP_UE(l,t,k)).*(H_eq_aux2(:,l)-HMean_AP_UE(l,t,k))';
                    end
                end
            end
            % Estimar el canal agregado
            [Hhat_agregated,B_agregated,C_agregated] = channelEstimates(H_eq,HMean_AP_UE,[],[],[],[],R_eq,[],nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,[]);
        else
            [Hhat_agregated,B_agregated,C_agregated] = channelEstimates(H_cascade,HMean_AP_UE,[],[],[],[],R_cascade,[],nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,[]);
            % H_eq = H_AP_UE;
            % R_eq = R_AP_UE;
        end
        % % Calcular SE
        [SE_P_MMSE, SE_MR_dist] = SE_uplink(Hhat_agregated,H_cascade,D,B_agregated,C_agregated,tau_c,tau_p,nbrOfRealizations,N_AP,K,L,p,R_cascade,pilotIndex);
        % 
        % 
        % Guardar resultados en la dimensión
        SE_PMMSE_DCC(:,n,s) = SE_P_MMSE;
        %sum(SE_P_MMSE)
        %SE_MR_DIST(:,n,s)  = SE_MR_dist;

        clear Hhat H_eq B C R_eq;
    end
end

save('results1')

%% Graficar resultados
figure; hold on; box on;
set(gca,'fontsize',16);

% P-MMSE 
aux1 = SE_PMMSE_DCC(:,:,1); % 0 RIS
aux2 = SE_PMMSE_DCC(:,:,2); % 5 RIS
aux3 = SE_PMMSE_DCC(:,:,3); % 10 RIS
aux4 = SE_PMMSE_DCC(:,:,4); % 20 RIS
aux5 = SE_PMMSE_DCC(:,:,5); % 50 RIS
aux6 = SE_PMMSE_DCC(:,:,6); % 100 RIS

plot(sort(aux1(:)), linspace(0,1,K*nbrOfSetups), 'k-', 'LineWidth', 2);
plot(sort(aux2(:)), linspace(0,1,K*nbrOfSetups), 'r-',  'LineWidth', 2);
plot(sort(aux3(:)), linspace(0,1,K*nbrOfSetups), 'g-', 'LineWidth', 2);
plot(sort(aux4(:)), linspace(0,1,K*nbrOfSetups), 'b-', 'LineWidth', 2);
plot(sort(aux5(:)), linspace(0,1,K*nbrOfSetups), 'm-',  'LineWidth', 2);
plot(sort(aux6(:)), linspace(0,1,K*nbrOfSetups), 'y-', 'LineWidth', 2);

% % Ejes y leyenda
xlabel('Spectral efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
ylabel('CDF', 'Interpreter', 'Latex');
legend({'P-MMSE 0 RIS', 'P-MMSE 5 RIS', 'P-MMSE 10 RIS', 'P-MMSE 20 RIS', 'P-MMSE 50 RIS','P-MMSE 75 RIS'}, 'Interpreter', 'Latex', 'Location', 'SouthEast');
xlim([0 25]);
