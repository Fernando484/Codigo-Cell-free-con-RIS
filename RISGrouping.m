function [H_cascade_grouped] = RISGrouping(H_unsummed, Ngroup, groupRIS_size, L,N_AP,S,nbrOfRealizations,K,N_RIS,N_V_RIS,N_H_RIS)

g_V = sqrt(groupRIS_size); % Elementos por fila en el grupo
g_H = sqrt(groupRIS_size); % Elementos por columna en el grupo

H_cascade_grouped = zeros(L*N_AP, Ngroup*S, nbrOfRealizations, K);

for k = 1:K
    for t = 1:nbrOfRealizations
        for m = 1:S
            % Extraer todos los elementos de ESTA RIS específica (m)
            idx_ris = (m-1)*N_RIS + 1 : m*N_RIS;

            for a = 1:L*N_AP
                % 1. Extraer el vector 1D de la cascada y pasarlo a grilla 2D (N_V_RIS x N_H_RIS)
                H_casc_2D = reshape(H_unsummed(a, idx_ris, t, k), [N_V_RIS, N_H_RIS]);

                % 2. Crear matriz temporal para el grupo
                H_group_temp = zeros(N_V_RIS/g_V, N_H_RIS/g_H);

                % 3. Recorrer la grilla y sumar los bloques 2x2
                for i = 1:(N_V_RIS/g_V)
                    for j = 1:(N_H_RIS/g_H)
                        bloque = H_casc_2D((i-1)*g_V+1 : i*g_V, (j-1)*g_H+1 : j*g_H);
                        H_group_temp(i, j) = sum(bloque, 'all'); % Suma los 4 elementos
                    end
                end

                % 4. Volver a aplanar (vector 1D) y asignar a la matriz final
                idx_group = (m-1)*Ngroup + 1 : m*Ngroup;
                H_cascade_grouped(a, idx_group, t, k) = reshape(H_group_temp, [1, Ngroup]);
            end
        end
    end
end