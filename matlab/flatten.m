function out = flatten(x, options)
%% Reduce dimensionality of input array
%  out = flatten(x, options)
%
% The size of the output depends on the input parameters. numel(out) will
% always equal numel(x), but the data is stored in different shapes using
% the reshape() function. The optional startDim and endDim arguments
% determine the target output size.
%
% startDim determines the first dimension to flatten, while endDim sets the
% last dimension to flatten. All dimensions in the range [startDim, endDim]
% inclusive are collapsed into a single dimension.
%
% Negative numbers may be used as shorthand for the last dimensions in the
% array, with -1 being the last dimension, -2 being the second to last
% dimension, and so on. 
%
% The default parameters are startDim=1, endDim=-1, which collapses all of
% the array's dimensions to create a column vector with length numel(x).
% 
% If startDim == endDim, only a single dimension is selected to be
% "flattened" and the output has the same size as the input.
% 
% Input
% -----
% x : matrix
%   array to flatten
%
% Options
% -------
% startDim : nonzero integer
%   The first dimension to flatten. Must be nonzero. Default is 1, the 
%   first dimension
% endDim : nonzero integer
%   The last dimension to flatten. Must be nonzero. Default is -1, the last 
%   dimension
% range : pair of nonzero integers
%   Shortcut for setting startDim and endDim, sets startDim to range(1) and 
%   endDim to range(end).
%
% Output
% ------
% out : matrix
%   flattened version of the input x
%
% Example
% -------
% >> a = zeros(1, 2, 3, 4, 5, 6, 7, 8);
% >> size(a)
% ans =
%   1   2   3   4   5   6   7   8
% >> size(flatten(a))
% ans =
%   40320   1
% >> size(flatten(a, 'endDim', -3))
% ans =
%   720   7   8
% >> size(flatten(a, 'startDim', 5))
% ans =
%   1   2   3   4   1680
% >> size(flatten(a, 'range', [4, 6]))
% ans =
%   1   2   3   120   7   8
%
%% Created 2022-07-15 by Samuel Adams-Tew
arguments
    x
    options.startDim = 1
    options.endDim = -1 % must be nonzero
    options.range
end

if isfield(options, 'range')
    startDim = options.range(1);
    endDim = options.range(end);
else
    startDim = options.startDim;
    endDim = options.endDim;
end

if startDim == 0
    error("startDim must be nonzero")
elseif startDim > ndims(x)
    startDim = -1;
end

if endDim == 0
    error("endDim must be nonzero")
elseif endDim > ndims(x)
    endDim = -1;
end

if startDim < 0
    startDim = ndims(x) + startDim + 1;
end
if endDim < 0
    endDim = ndims(x) + endDim + 1;
end

if startDim > endDim
    error("If startDim and endDim are different, startDim must come before endDim")
elseif startDim == endDim
    % Output is same shape as input, return early
    out = x; return
end

% Determine the dimensions that will be preserved
keep1 = 1:(startDim-1);
if ~isempty(keep1); keep1 = num2cell(size(x, keep1)); end
keep2 = (endDim+1):ndims(x); 
if ~isempty(keep2); keep2 = num2cell(size(x, keep2)); end

% Apply reshape, using the determined preserved dimensions
if isempty(keep1) && isempty(keep2)
    out = x(:);
elseif isempty(keep1) && ~isempty(keep2)
    out = reshape(x, [], keep2{:});
elseif ~isempty(keep1) && isempty(keep2)
    out = reshape(x, keep1{:}, []);
else
    out = reshape(x, keep1{:}, [], keep2{:});
end

end