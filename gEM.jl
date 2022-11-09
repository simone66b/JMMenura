using LinearAlgebra, ExponentialUtilities

P0 = [0.6856935  -0.0424340    1.2743154   0.5162707; 
-0.0424340   0.1899794   -0.3296564  -0.1216394; 
1.2743154  -0.3296564    3.1162779   1.2956094; 
0.5162707  -0.1216394    1.2956094   0.5810063];

G = ([0 1 0 0; -1 0 0 0; 0 0 0 0; 0 0 0 0], 
[0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0], 
[0 0 0 1; 0 0 0 0; 0 0 0 0; -1 0 0 0], 
[0 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 0],
[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 -1 0 0], 
[0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0]); ## basis of skew-symmetric 4x4-matrices

a = randn(6); b = randn(6); ## these values can be set arbitrarily, even constant
A = sum(a.*G); B = sum(b.*G); ## random skew-symmetric matrices

dt = 0.1;
t = dt:dt:1; 
N = length(t);

M = 2;	## number of different paths
dW = randn(N,M)*sqrt(dt); ## Brownian motion

P = zeros(4,4,N+1,M);

for m = 1:M
	P[:,:,1,m] = P0;
	for n = 1:N
		AA = t[n]*A; BB = t[n]*B; ## time-dependent drift and diffusion coefficients
		Omega = AA*dt + BB*dW[n,m]; ## Euler-Maruyama
		expO = exponential!(Omega);
		P[:,:,n+1,m] = expO*P[:,:,n,m]*expO';
	end
end