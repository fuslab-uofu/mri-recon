function w0 = ndseparablewindow(windowfun, sz, param)
%% Window function with arbitrary dimensions.
%  w0 = ndseparablewindow(windowfun, sz, param)
% 
% Generates a multidimensional window as the outer product of an arbitrary
% number of 1-D windows.
%
% Recommended to use with window functions in the MATLAB Signal Processing
% Toolbox.
%
% See also TUKEYWIN, GAUSSWIN, CHEBWIN, KAISER, FERMIWIN.
%
% Input
% -----
% windowfun : function handle, (sz: int, param: scalar) -> column vector
%   MATLAB 1D window function (such as @tukeywin); in the function
%   signature specification, `sz` is the number of points, `param` is some 
%   characteristic of the window (usually related to bandwidth), and the
%   output is a 1D window column vector.
% sz : vector of positive integers
%   The size of each dimension of the filter
% param : 
%   The parameter to pass for each dimension (in same order as sz)
% 
% Output
% ------
% w0 : matrix
%   An window function with size defined by input `sz`, `size(w0) == sz`.
%
% Example
% -------
% >> ndseparablewindow(@tukeywin, [4, 8], [1, 1])
% ans =
%   0   0       0       0       0       0       0       0
%   0   0.1412  0.4584  0.7129  0.7129  0.4584  0.1412  0
%   0   0.1412  0.4584  0.7129  0.7129  0.4584  0.1412  0
%   0   0       0       0       0       0       0       0
%
%% Created 2022-04-05 Samuel Adams

ndims = length(sz);
axes = cell(1, ndims);
grids = cell(1, ndims);
for dim = 1:length(sz)
    if sz(dim) == 1
        % singleton dimensions have value of 1 everywhere
        axes{dim} = 1;
    else
        axes{dim} = windowfun(sz(dim), param(dim));
    end
end

[grids{:}] = ndgrid(axes{:});

gridmat = cat(ndims + 1, grids{:});
w0 = prod(gridmat, ndims + 1);

end