function im = ifftc(ks, dim, n)
%% Centered inverse fast Fourier transform
%  ks = ifftc(im, dim, n)
% 
% Performs a centered (n-point) inverse fast fourier transform along the 
% specified dimension(s). Equivalent to performing
%       fftshift(ifft(ifftshift(ks, dim), n, dim), dim)
% along each of the given dimensions.
%
% Input
% -----
% ks : matrix
%   K-space data with arbitrary dimensions.
% dim : vector of positive integers, optional
%   Dimensions along which to apply inverse Fourier transform. If not 
%   specified, applied to all dimensions.
% n : vector of positive integers, optional
%   Number of points to use in each Fourier transform specified by dim.
%   If not specified, n = size(im, dim) (Fourier transform has same
%   number of points as input). See IFFT documentation for more 
%   information.
%
% Output
% ------
% im : numeric matrix
%   Inverse fourier transformed (image-space) data. If n is not 
%   specified, has same size as the input. Otherwise, any dimensions 
%   included in dim have the number of points specified in n.
%
% See also IFFT.
%
%% Created 2023-10-13 Samuel Adams-Tew

if ~exist('dim', 'var')
    dim = 1:ndims(ks);
end

if exist('n', 'var')
    if length(n) == 1
        n = n*ones(size(dim));
    end
else 
    n = size(ks, dim);
end

im = ks;
for d = 1:length(dim)
    im = fftshift(ifft(ifftshift(im, dim(d)), n(d), dim(d)), dim(d));
end

end