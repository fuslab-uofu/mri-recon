function cropIm = crop_to_size(im, szOut)
%% Crop image data
%  cropIm = crop_to_size(im, szOut)
% 
% Crops the input image to the given size using symmetric cropping. If an 
% odd number of entries are to be removed along a given dimension, the 
% extra entry is removed from the front of that dimension. Can be used to
% exactly undo the action of pad_to_size().
%
% Input
% -----
% im : matrix
%   Image data with arbitrary dimensions
% szOut : vector of positive integers
%   The desired output size. If length(szOut) < ndims(im), dimensions with 
%   no size specified in szOut will be left as-is. if 
%   szOut(dim) > size(im, dim), then the size for that dimension will also 
%   be left as-is.
%
% Output
% ------
% cropIm : matrix 
%   Cropped image with same number of dimensions as `im` and the size in
%   each dimension dim <= length(szOut) matching the requested size.
%
% Example
% -------
% >> cropIm = crop_to_size(zeros([8, 8, 8, 8]), [3, 4, 5])
% >> size(cropIm)
% ans =
%   3   4   5   8
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
    if diff < 0
        % If output is shorter than input, create ranges that induce
        % truncation
        truncS = ceil(-diff/2);
        truncE = -diff - truncS;
        
        sourceRng{d} = (1 + truncS):(szIn(d) - truncE);
        destRng{d} = 1:szOut(d);
    else
        if diff > 0
            warning('Requested size for dimension %d is larger than input. Maintaining input size for dimension %d.', d, d)
            szOut(d) = szIn(d);
        end
        % Create ranges that maintain size
        sourceRng{d} = 1:szIn(d);
        destRng{d} = 1:szOut(d);
    end
end

% copy values into new image using ranges found above
cropIm = zeros(szOut, 'like', im);
cropIm(destRng{:}) = im(sourceRng{:});

end