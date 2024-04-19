function imZFI = zfi_to_size(im, szOut)
%% Zero filled interpolation
%  imZFI = zfi_to_size(im, szOut)
%
% Interpolates the image to the specified size using zero padding in the
% spatial frequency domain.
%
% Intepolation is achieved by applying the Fourier transform along each 
% dimension specified in szOut. If szOut(dim) > size(im, dim), 
% fft(im, [], dim) is padded symmetrically with zeros to have size
% szOut(dim). szOut(dim) > size(im, dim), fft(im, [], dim) is symmetrically
% truncated to have size szOout(dim).
%
% Input
% -----
% im : matrix 
%   Image data with arbitrary number of dimensions
% szOut : vector of positive integers
%   The desized output size. If length(szOut) < ndims(im), dimensions with 
%   no size specified in szOut will be left as-is.
%
% Output
% ------
% imZFI : complex matrix
%   Complex valued matrix with zero-filled interpolation to specified size.
%
%% Created 2022-07-15 Samuel Adams-Tew

% Get number of output dimensions
ndim = max(ndims(im), length(szOut));
% Get current size
szIn = size(im, 1:ndim);
% Determine which dimensions are being ZFI'd
zfiDims = find(szIn < szOut);

% Apply ZFI process
% FT -> zero pad -> inverse
ks = fftc(im, zfiDims);
ksZFI = pad_to_size(ks, szOut);
imZFI = ifftc(ksZFI, zfiDims);

end