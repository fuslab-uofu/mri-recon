function imPad = pad_to_size(im, szOut, padVal)
%% padIm = pad_to_size(im, szOut)
% Pads the input image to a given size using symmetric padding. If an odd
% number of entries are to be added along a given dimension, the extra
% entry is added at the front of that dimension.
%
% Input
% -----
% im : matrix 
%   Image data with arbitrary number of dimensions
% szOut : vector of positive integers
%   The desized output size. If length(szOut) < ndims(im), dimensions with 
%   no size specified in szOut will be left as-is.
% padVal : scalar, optional
%   Value to use for padding the image to the requested size. Default is 0.
%
% Output
% ------
% imPad : matrix
%   Image padded to the requested size with im at the center and padVal
%   everywhere else.
%
% Examples
% --------
% >> a = [1, 2, 3; 4, 5, 6]
% a =
%      1     2     3
%      4     5     6
% >> pad_to_size(a, [5, 5])
% ans =
%      0     0     0     0     0
%      0     0     0     0     0
%      0     1     2     3     0
%      0     4     5     6     0
%      0     0     0     0     0
% >> pad_to_size(a, [4, 6], -10)
% ans =
%    -10   -10   -10   -10   -10   -10
%    -10   -10     1     2     3   -10
%    -10   -10     4     5     6   -10
%    -10   -10   -10   -10   -10   -10
% >> size(pad_to_size(a, 8))
% ans =
%      8     3
%
%% 2023-10-13 Samuel Adams-Tew

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
        % If output is shorter than input, create ranges that maintain size
        % and warn the user
        warning('Requested size for dimension %d is smaller than input. Maintaining input size for dimension %d.', d, d)
        sourceRng{d} = 1:szIn(d);
        destRng{d} = 1:szOut(d);
    else
        % Create ranges that maintain size
        sourceRng{d} = 1:szIn(d);
        destRng{d} = 1:szOut(d);
    end
end

% Create array of padVal at the new size
if exist('padVal', 'var')
    imPad = padVal * ones(szOut, 'like', im);
else
    imPad = zeros(szOut, 'like', im);
end

% Copy values into padded image
imPad(destRng{:}) = im(sourceRng{:});

end