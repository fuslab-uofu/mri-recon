function imZFI = zfi_to_size(im, szOut)
%% imZFI = zfi_to_size(im, szOut)
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

% Get current size
szIn = size(im);

sourceRng = cell(1, ndims(im));
destRng = cell(1, ndims(im));
% Calculate padding amounts along each dimension
for d = 1:ndims(im)
    if length(szOut) < d
        szOut(d) = szIn(d);
    end
    % Get difference between input and output size
    diff = szOut(d) - szIn(d);
    if diff > 0
        % If output is longer than input, create ranges that induce padding
        padS = ceil(diff/2);
        padE = diff - padS;

        sourceRng{d} = 1:szIn(d);
        destRng{d} = (1 + padS):(szOut(d) - padE);
    elseif diff < 0
        % If output is shorter than input, create ranges that induce
        % truncation
        truncS = ceil(-diff/2);
        truncE = -diff - truncS;

        sourceRng{d} = (1 + truncS):(szIn(d) - truncE);
        destRng{d} = 1:szOut(d);
    else
        % Create ranges that maintain size
        sourceRng{d} = 1:szIn(d);
        destRng{d} = 1:szOut(d);
    end
end

% Perform FT along axes that will change size
ksIn = im;
for d = 1:ndims(im)
    if szOut(d) ~= szIn(d)
        ksIn = fftshift(fft(ifftshift(ksIn, d), [], d), d);
    end
end

% copy values into new ZFI k-space using ranges found above
ksOut = zeros(szOut, 'like', ksIn);
ksOut(destRng{:}) = ksIn(sourceRng{:});

% Perform IFT along axes that change size
imZFI = ksOut;
for d = 1:ndims(im)
    if szOut(d) ~= szIn(d)
        imZFI = fftshift(ifft(ifftshift(imZFI, d), [], d), d);
    end
end

end