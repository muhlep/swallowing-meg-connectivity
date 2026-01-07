function [cleanedData, artLog] = ft_CleanTrials(data, varargin)
% ft_CleanTrials
% Deterministic, batch-friendly trial rejection for MEG/EEG raw data.
% Rejects trials based on:
%   1) absolute amplitude (max |x|)
%   2) robust z-score (max |(x-median)/MAD|)
%   3) temporal gradient (max |diff(x)|)
%   4) presence of NaNs/Infs
%
% QC gate:
%   If fewer than MinTrials remain, returns an EMPTY dataset (trial/time cleared).
%
% INPUT
%   data : FieldTrip raw structure (data.trial{t} = [nChan x nTime])
%
% NAME–VALUE PAIRS
%   'AmpThreshold'    : absolute amplitude threshold (default: 6e-12)
%   'ZThreshold'      : robust z-score threshold using MAD (default: 10)
%   'GradThreshold'   : temporal gradient threshold (default: 6e-12)
%   'MinTrials'       : minimum trials required (default: 5)
%   'Verbose'         : true/false (default: true)
%
% OUTPUTS
%   cleanedData : same structure as input with bad trials removed
%   artLog      : struct with rejection summary + per-trial metrics
%
% Notes:
% - Uses robust MAD-based z-score (more stable than mean/std in presence of outliers).
% - Does NOT rebuild sampleinfo. Lets FieldTrip handle it via ft_checkdata.
% - Keeps data.label/time, removes trial/time entries consistently.

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('AmpThreshold',  6e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ZThreshold',    10.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('GradThreshold', 6e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinTrials',     5,     @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Verbose',       true,  @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

assert(isfield(data,'trial') && iscell(data.trial) && ~isempty(data.trial), ...
    'Input must be FieldTrip raw data with non-empty data.trial.');

nTrials = numel(data.trial);

% -------------------- preallocate metrics --------------------
maxAmp      = nan(1, nTrials);
maxZrobust  = nan(1, nTrials);
maxGrad     = nan(1, nTrials);
hasBadVals  = false(1, nTrials);

% -------------------- compute metrics per trial --------------------
for t = 1:nTrials
    x = data.trial{t};                 % [nChan x nTime]
    xFlat = x(:);

    % NaN/Inf check
    hasBadVals(t) = any(~isfinite(xFlat));

    % Max absolute amplitude
    maxAmp(t) = max(abs(xFlat), [], 'omitnan');

    % Robust z-score: (x - median) / (1.4826*MAD)
    medx = median(xFlat, 'omitnan');
    madx = median(abs(xFlat - medx), 'omitnan');
    if madx < eps
        maxZrobust(t) = 0;
    else
        z = (xFlat - medx) / (1.4826 * madx);
        maxZrobust(t) = max(abs(z), [], 'omitnan');
    end

    % Temporal gradient max (across channels & time)
    dx = diff(x, 1, 2);
    maxGrad(t) = max(abs(dx), [], 'all', 'omitnan');
end

% -------------------- decide bad trials --------------------
badIdx = (maxAmp     > opt.AmpThreshold) | ...
         (maxZrobust > opt.ZThreshold)   | ...
         (maxGrad    > opt.GradThreshold)| ...
         hasBadVals;

removedTrials = find(badIdx);
nRemoved = numel(removedTrials);
nKept = nTrials - nRemoved;

% -------------------- logging --------------------
artLog = struct();
artLog.nOrig   = nTrials;
artLog.nRejected = nRemoved;
artLog.removed = removedTrials(:)';
artLog.criteria = struct( ...
    'AmpThreshold',  opt.AmpThreshold, ...
    'ZThreshold',    opt.ZThreshold, ...
    'GradThreshold', opt.GradThreshold, ...
    'MinTrials',     opt.MinTrials);
artLog.metrics = struct( ...
    'maxAmp',     maxAmp, ...
    'maxZrobust', maxZrobust, ...
    'maxGrad',    maxGrad, ...
    'hasBadVals', hasBadVals);

if opt.Verbose
    fprintf('\n[ft_CleanTrials] Removed %d/%d trials (%.1f%%). Kept: %d\n', ...
        nRemoved, nTrials, 100*nRemoved/max(nTrials,1), nKept);
end

% -------------------- apply rejection --------------------
cleanedData = data;
cleanedData.trial(badIdx) = [];
if isfield(cleanedData,'time') && numel(cleanedData.time)==nTrials
    cleanedData.time(badIdx) = [];
end

% Remove sampleinfo rows if present
if isfield(cleanedData,'sampleinfo') && size(cleanedData.sampleinfo,1)==nTrials
    cleanedData.sampleinfo(badIdx,:) = [];
end

% -------------------- QC gate --------------------
if nKept < opt.MinTrials
    warning('[ft_CleanTrials][QC] Fewer than %d trials remain (%d/%d). Returning empty dataset.', ...
        opt.MinTrials, nKept, nTrials);
    cleanedData.trial = {};
    if isfield(cleanedData,'time'), cleanedData.time = {}; end
    if isfield(cleanedData,'sampleinfo'), cleanedData.sampleinfo = []; end
    return;
end

% -------------------- make FT-consistent --------------------
% This will fix minor inconsistencies (including sampleinfo/time alignment).
try
    cleanedData = ft_checkdata(cleanedData, 'datatype','raw', 'hassampleinfo','yes');
catch
    % If ft_checkdata not available (should be), just return as is.
end

end