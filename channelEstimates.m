function [Hhat,B,C,Hhat_cascade,B_cascade,C_cascade] = channelEstimates(H_AP_UE,HMean_AP_UE,HMean_AP_RIS,HMean_RIS_UE,H_AP_RIS,H_RIS_UE,R_AP_UE,R_cascade,nbrOfRealizations,L,K,N_AP,tau_p,pilotIndex,p,risAssignment,S,Ngroup)
% Esta función estima el canal directo entre los UEs y los APs y además
% estima los canales en cascada entre UE-RIS-AP usando el estimador MMSE
% Argumentos de entrada:
% H_AP_UE                   = Canal Rician generado entre el AP y el UE de
%                             dimensiones (L*N_AP,nbrOfRealizations,K)
% HMean_AP_UE               = Componente de LoS del canal AP-UE
% HMean_AP_RIS              = Componente de LoS del canal AP-RIS
% HMean_RIS_UE              = Componente de LoS del canal RIS-UE
% H_AP_RIS                  = Canal Rician generado entre el AP y la RIS de
%                             dimensiones (L*N_AP,S*N_RIS,nbrOfRealizations)
% H_RIS_UE                  = Canal Rician generado entre la RIS y el UE de
%                             dimensiones (S*N_RIS,nbrOfRealizations, K)
% R_AP_UE                   = Matriz de correlación espacial entre el AP y el UE
%                             de dimensiones (L*N_AP,L*N_AP,L,K)
% R_cascade                 = Matriz de correlación espacial del canal en
%                             cascada AP-RIS-UE de dimensiones (N_AP,N_AP,L,K)
% nbrOfRealizations         = Número de realizaciones del canal  
% L                         = Número de APs
% K                         = Número de usuarios
% N_AP                      = Número de antenas de cada AP
% N_RIS                     = Número de elementos de la RIS
% tau_p                     = Número de símbolos de un bloque de coherencia
%                             dedicados a una secuencia de pilotos
% pilotIndex                = Matriz que indica la secuencia de pilotos
%                             asignada a cada UE
% p                         = Potencia de transmisión en mW
% risAssignment             = Cell array que indica que UEs están asignados
%                             a cada RIS
% S                         = Número de RIS
% Ngroup                    = Número de agrupaciones de elementos de una
%                             RIS
% Argumentos de salida:
% Hhat                      = Matriz de los canales directos estimados
% Hhat_cascade              = Matriz de los canales en cascada estimados
% B                         = Correlación de estimación del canal directo
% C                         = Correlación del error de estimación del canal
%                             directo
% B_cascade                 = Correlación de estimación del canal en
% cascada
% C_cascade                 = Correlación del error de estimación del canal
%                             en cascada

%% Estimación del canal

% Matriz identidad de tamaño N_APxN_AP
eyeN_AP = eye(N_AP);

% Ruido normalizado para la estimación de piloto AP-UE
Np = sqrt(0.5)*(randn(N_AP,nbrOfRealizations,L,tau_p) + 1i*randn(N_AP,nbrOfRealizations,L,tau_p));

% Preparar para almacenar los resultados AP-UE
Hhat = zeros(L*N_AP,nbrOfRealizations,K);

% Preparar para almacenar los resultados AP-RIS-UE
if (nargout>3)
    Hhat_cascade = zeros(L*N_AP,S*Ngroup,nbrOfRealizations,K);
    H_eq_tl = zeros(L*N_AP,nbrOfRealizations);
end

% Reservar matriz de correlación de estimación si se requiere
if nargout>1
    B = zeros(size(R_AP_UE));
end

% Reservar matriz de error de estimación si se requiere
if nargout>2
    C = zeros(size(R_AP_UE));
end

% Reservar matriz de correlación de estimación si se requiere
if nargout>4
    B_cascade = zeros(N_AP,N_AP,L,K);
end

% Reservar matriz de error de estimación si se requiere
if nargout>4
    C_cascade = zeros(N_AP,N_AP,L,K);
end

% Para cada AP
for l = 1:L
    
    % Para cada piloto
    for t = 1:tau_p
        
        % Señal recibida procesada
        yp = sqrt(p)*tau_p*sum(H_AP_UE((l-1)*N_AP+1:l*N_AP,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
        yMean = zeros(N_AP, nbrOfRealizations);  % Inicializar la media esperada de y
        if (nargout>3)
            yMean_cascade = zeros(N_AP, nbrOfRealizations);
        end

        % Matriz a invertir en el estimador MMSE
        PsiInv = (p*tau_p*sum(R_AP_UE(:,:,l,t==pilotIndex),4) + eyeN_AP);
        
        % Para cada usuario que usa el piloto t
        for k = find(t==pilotIndex)'
            % Calcular estimación del canal directo
            RPsi = R_AP_UE(:,:,l,k) / PsiInv;
            yMean = yMean + sqrt(p)*tau_p*HMean_AP_UE((l-1)*N_AP+1:l*N_AP,:,k);
            Hhat((l-1)*N_AP+1:l*N_AP,:,k) = sqrt(p)*RPsi*(yp - yMean) + HMean_AP_UE((l-1)*N_AP+1:l*N_AP,:,k);

            % Calcular la estimación del canal en cascada hacia las RIS que
            % tiene asignado cada usuario
            if (S> 0 && ismember(k, [risAssignment{:}]) && nargout>3)                            % Buscar si hay alguna ris asignada a un usuario
                ris_k = find([risAssignment{:}]==k);                        % Ver que ris hay asignadas
                for i = 1:length(ris_k)                                     % Recorrer las ris asignadas
                    s = ris_k(i);
                    for j = 1:Ngroup                                         % Recorrer cada elemento individual de la RIS
                        %R_cascade = R_AP_RIS1(1:N_AP,1:N_AP,l,s)*R_AP_RIS2(j,j,l,s)*R_RIS_UE(j,j,s,k);
                        PsiInv_cascade = (p*tau_p*sum(R_cascade(:,:,l,t==pilotIndex),4) + eyeN_AP);
                        RPsi_cascade = R_cascade(:,:,l,k)/PsiInv_cascade;
                        for r = 1:nbrOfRealizations
                            % Canal equivalente suma del canal directo más
                            % el canal en cascada para un elemento de la
                            % RIS encendido con fase 0
                            H_eq_tl((l-1)*N_AP+1:l*N_AP,r) = (H_AP_RIS((l-1)*N_AP+1:l*N_AP,j + (Ngroup * (s - 1)),r) * squeeze(H_RIS_UE(j + (Ngroup * (s - 1)), r, k)));
                            yMean_cascade(:,r) = yMean_cascade(:,r) + sqrt(p)*tau_p*HMean_AP_RIS((l-1)*N_AP+1:l*N_AP,j + (Ngroup * (s - 1)),r)*HMean_RIS_UE(j + (Ngroup * (s - 1)),r,k);
                        end
                        %yMean_cascade = yMean_cascade + sqrt(p)*tau_p*HMean_AP_RIS((l-1)*N_AP+1:l*N_AP,:,k)*HMean_RIS_UE((l-1)*N_AP+1:l*N_AP,:,k);
                        yp_cascade = sqrt(p)*tau_p*H_eq_tl((l-1)*N_AP+1:l*N_AP,:) + sqrt(tau_p)*Np(:,:,l,t) - yp; % Señal recibida por el canal en cascada de los elementos activos de la ris
                        Hhat_cascade((l-1)*N_AP+1:l*N_AP,j + (Ngroup * (s - 1)),:,k) = sqrt(p)*RPsi_cascade*(yp_cascade-yMean_cascade) + HMean_AP_RIS((l-1)*N_AP+1:l*N_AP,j + (Ngroup * (s - 1)),r)*HMean_RIS_UE(j + (Ngroup * (s - 1)),r,k);
                        % Correlación de estimación canal en cascada
                        if nargout>4
                            B_cascade(:,:,l,k) = p*tau_p*RPsi_cascade*R_cascade(:,:,l,k);
                        end
                        
                        % Correlación del error de estimación canal en
                        % cascada
                        if nargout>4
                            C_cascade(:,:,l,k) = R_cascade(:,:,l,k) - B_cascade(:,:,l,k);
                        end
                    end
                end
            end
            
            % Correlación de estimación canal directo
            if nargout>1
                B(:,:,l,k) = p*tau_p*RPsi*R_AP_UE(:,:,l,k);
            end
            
            % Correlación del error de estimación directo
            if nargout>2
                C(:,:,l,k) = R_AP_UE(:,:,l,k) - B(:,:,l,k);
            end
            
        end
        
    end
    
end

