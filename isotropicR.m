function R = isotropicR(lambda, N_H, N_V, d) 
% lambda es la longitud de onda
% N_H es el número de elementos de cada fila
% N_V es el número de elementos de cada columna
% d es el tamaño de cada elemento cuadrado de la RIS dH x dV

%Generate a grid for the elements
gridPoints_H = (0:N_H-1)*d;
gridPoints_V = (0:N_V-1)*d;


[X,Y] = meshgrid(gridPoints_H,gridPoints_V);

locations = X(:)+1i*Y(:);


%Total number of elements
N = length(locations);


%Compute the spatial correlation matrix
R = zeros(N,N);

for m = 1:N
    for l = 1:N
        
        R(m,l) = sinc(2*abs(locations(m)-locations(l))/lambda);
        
    end
end