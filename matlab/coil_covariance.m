function [covmat, invcov] = coil_covariance(noise, dims)
%% Coil noise covariance matrix
%  [covmat, invcov] = coil_covariance(noise, dims)
%
% Computes the coil covariance based on noise data.
%
% Input
% -----
% noise : complex matrix
%   Coil noise data.
% dims : cell array of string
%   Array with the dimension labels for each dimension in `noise`. Must
%   include 'CHA' label, signifying channels / coils.
%
% Output
% ------
% covmat : square matrix
%   Covariance matrix derived from coil noise, where each entry is the
%   correlation between pairs of coil elements (diagonal is
%   self-correlation). Square matrix with nCha rows / columns.
% invcov : square matrix
%   Inverse of `covmat`.
%
%% Created 2024-04-09 Samuel Adams-Tew

% Argument validation
if length(dims) ~= ndims(noise)
    error("dim must have exactly one entry for each dimension in noise")
end

% Permute noise data to have channels as first dimension
order = get_permutation(dims, {'CHA'});
noise = permute(noise, order);

% Flatten noise data to 2D, with channels along column dimension
noise = transpose(flatten(noise, 'startDim', 2));

% Compute covariance and inverse
covmat = cov(noise); invcov = inv(covmat);

end