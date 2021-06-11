function x = ifft2c(x)
%
% res = ifft2c(x)
% 
% Orthonormal inverse 2D FFT on the first 2 dimensions of the input array.
% 
% (c) Kangrong Zhu,      Stanford University    2011

sz = size(x);
n_2d = prod( sz(3 : end) ); % # of 2D arrays to perform 2D iFFT on.

x = reshape(x, [sz(1:2), n_2d]);
x = sqrt(prod(sz(1:2))) * fftshift(ifft2(ifftshift(x)));
x = reshape(x, sz);
