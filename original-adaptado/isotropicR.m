function R = isotropicR(lambda, N_HV, d) 
% lambda es la longitud de onda
% N_HV es el número de elementos de cada fila y cada columna
% d es el tamaño de cada elemento cuadrado de la RIS dH x dV

%Generate a grid for the elements
gridPoints = (0:N_HV-1)*d;

[X,Y] = meshgrid(gridPoints,gridPoints);

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