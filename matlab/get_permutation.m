function [permuteOrder, permuteLabels] = get_permutation(current, new, options)
%% [permuteOrder, permuteLabels] = get_permutation(current, new, options)
% Given dimension label cell arrays, determines the dimension order that
% produces the desired ordering
%
% `permuteOrder` contains the dimension numbers that match the ordering in
% new. If there are entries in current that are not in new, these are put
% at the end, such that all dimensions specified in new are in the 
% specified order at the beginning, and all unspecified dimensions follow
% at the end.
%
% Input
% -----
% current : cell array of strings
%   Dimension labels for the current dimension order
% new : cell array of strings
%   Desired dimension ordering
%
% Options
% -------
% ignoreCase : logical
%   Sets case sensitive compare flag for comparing dimension labels.
%   Default is true.
% insertMissing : logical
%   Whether to insert dimensions if entries in new are not included in 
%   `current`. If true, singleton dimensions will be inserted for
%   each entry in new that is not in current. Otherwise, these dimensions
%   are omitted. Default is true.
%
% Output
% ------
% permuteOrder : vector of positive integers
%   The dimension order to pass to the permute() function to order the 
%   array as specified.
% permuteLabels : cell array of strings
%   New dimension labels after the array is permuted as specified.
%
% Example
% -------
% >> dims = {'COL' 'CHA' 'LIN' 'ACQ'};
% >> get_permutation(dims, {'lin', 'col', 'cha', 'acq'})
% ans =
%   3   1   2   4
% >> get_permutation(dims, {'cha', 'acq', 'ech'})
% ans =
%   2   4   5   1   3
% >> get_permutation(dims, {'cha', 'acq', 'ech'}, 'insertMissing', false)
% ans =
%   2   4   1   3
%
%% Created 2023-01-13 Samuel Adams-Tew
arguments
    current
    new
    options.ignoreCase = true
    options.insertMissing = true
end

permuteOrder = [];
permuteLabels = {};
unusedDims = 1:length(current);
missingDims = [];

for d = 1:length(new)
    idx = find(matches(current, new{d}, 'IgnoreCase', options.ignoreCase));
    if isempty(idx)
        missingDims = [missingDims d];
    else
        permuteOrder = [permuteOrder idx];
        permuteLabels = [permuteLabels current(idx)];
        unusedDims(unusedDims == idx) = [];
    end
end

permuteOrder = [permuteOrder unusedDims];
permuteLabels = [permuteLabels current(unusedDims)];

insertNum = length(current) + 1;
if options.insertMissing && ~isempty(missingDims)
    insertNum = length(current) + 1;
    for d = 1:length(missingDims)
        missing = missingDims(d);
        if missing == 1
            permuteOrder = [insertNum permuteOrder];
            permuteLabels = [new(missing) permuteLabels];
        else
            permuteOrder = [permuteOrder(1:(missing-1)), insertNum, permuteOrder(missing:end)];
            permuteLabels = [permuteLabels(1:(missing-1)), new(missing), permuteLabels(missing:end)];
        end
        insertNum = insertNum + 1;
    end
end

end