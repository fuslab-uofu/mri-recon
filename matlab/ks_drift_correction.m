function kscorr = ks_drift_correction(ks, dims, options)
%% K-space based B0 field drift correction
%  kscorr = ks_drift_correction(ks, dim, nBaseline)
%
% Corrects for phase change due to B0 field drift in repeated Cartesian 
% k-space measurements, based on Ref. [1].
%
% Input
% -----
% ks : complex matrix
%   Coil noise data.
% dims : cell array of string
%   Array with the dimension labels for each dimension in `ks`.
%
% Options
% -------
% nBaseline : nonnegative integer
%   Number of baseline acquisitions in time. Default is 5.
% TimeDim : string
%   Label of dimension that contains k-space samples separated in time.
%   Default is 'REP'
% AverageDims: cell array of string
%   Labels of dimensions along which to average phase measurements. Default 
%   is {'COL', 'PAR', 'CHA'}.
% 
% Output
% ------
% kscorr : complex matrix
%   k-space data with B0 phase drift correction applied
%
% References
% ----------
% [1] Parker et al., "A k-space-based method to measure and correct for 
%   temporal B0 field variations in MR temperature imaging," doi: 
%   10.1002/mrm.29275
%
%% Updated 2024-04-16. Authors: Samuel Adams-Tew, Dennis L. Parker
arguments
    ks
    dims
    options.nBaseline = 5
    options.TimeDim = 'REP'
    options.AverageDims = {'COL', 'PAR', 'CHA'}
end

% Get index of time dimension
timeDim = find(contains(dims, 'REP', 'IgnoreCase', true));

% Create a cell array that indexes all entries along every dimension
idx = repmat({':'}, [1, ndims(ks)]);
% Replace timeDim entry with range through nBaseline
idx{timeDim} = 1:options.nBaseline;

% Mean value over time using indices defined above
ksmean = mean(ks(idx{:}), timeDim);
% Complex difference between actual samples and baseline average
ksdiff = ks .* conj(ksmean);

% Average differences over requested dimensions, compute phase angle
avgDims = find(contains(dims, options.AverageDims, 'IgnoreCase', true));
ksangcor = angle(mean(ksdiff, avgDims));

% Apply phase correction to original data
kscorr = ks .* exp(-1i*ksangcor);

end