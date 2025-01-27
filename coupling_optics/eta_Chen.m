function [etaF,FA] = eta_Chen(EA,w0,lambda,f,D)
%ETA_CHEN Computes coupling efficiency according to formula found on paper
%"Experimental demonstration of single-mode fiber coupling over relatively
%strong turbulence with adaptive optics" - Mo Chen,...

% Crear mesh conveniente con el tamaño tel.D para las
% variables X e Y para el back propagated mode (fórmula Chen).
sizeEA = size(EA,1);
x = (1:1:sizeEA) - 0.5*sizeEA;
x = D*x/sizeEA;
y = x;

[X Y]=meshgrid(x);
a = w0/lambda/f;
FA = sqrt(2*pi)*a.*exp(-(X.^2+Y.^2)*pi^2*a^2); % Back propagated mode
P = EA.*conj(FA);

% Integración numérica
% B = (trapz(X,trapz(Y, abs(P), 1))).^2;
% B = (trapz(X,trapz(Y, P, 1))).^2;
B = abs((trapz(x,trapz(y, P, 1))).^2);
C = trapz(x,trapz(y, abs(EA).^2, 1));
% F = trapz(X,trapz(Y, abs(FA).^2, 1));

etaF = B/C;

end

