function ks = fftc(im, dim, n)
%% Centered fast Fourier transform
%  ks = fftc(im, dim, n)
%
% Performs a centered (n-point) fast fourier transform along the specified
% dimension(s). Equivalent to performing
%       fftshift(fft(ifftshift(ks, dim), n, dim), dim)
% along each of the given dimensions.
%
% Input
% -----
% im : matrix
%   Image space data with arbitrary dimensions.
% dim : vector of positive integers, optional
%   Dimensions along which to apply Fourier transform. If not specified,
%   applied to all dimensions.
% n : vector of positive integers, optional
%   Number of points to use in each Fourier transform specified by dim.
%   If not specified, n = size(im, dim) (Fourier transform has same
%   number of points as input). See FFT documentation for more information.
%
% Output
% ------
% ks : numeric matrix
%   Fourier transformed data. If n is not specified, has same size as the
%   input. Otherwise, any dimensions included in dim have the number of
%   points specified in n.
%
% See also FFT.
%
%% Created 2023-10-13 Samuel Adams-Tew

if ~exist('dim', 'var')
    dim = 1:ndims(im);
end

if exist('n', 'var')
    if length(n) == 1
        n = n*ones(size(dim));
    end
else 
    n = size(im, dim);
end

ks = im;
for d = 1:length(dim)
    ks = fftshift(fft(ifftshift(ks, dim(d)), n(d), dim(d)), dim(d));
end

end