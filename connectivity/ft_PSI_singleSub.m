function psi_outfile = ft_PSI_singleSub(subj, cond, freqLabel, freqRange, varargin)
% ft_PSI_singleSub
% Compute subject-level PSI (Phase Slope Index) for one subject/condition/frequency band.
% Expects ROI time series saved by ft_Beamformer_extractROI:
%   roi_timeseries_<cond>_<freqLabel>.mat  with roi_data [nROI x nTrials x nSamples]
%
% INPUTS
%   subj       : subject ID string
%   cond       : condition label, e.g., 'EMG2', 'EMG3'
%   freqLabel  : short band label used in filenames, e.g., 'alpha','beta','lowGamma','highGamma','theta'
%   freqRange  : [fmin fmax] in Hz (defines PSI frequency support)
%
% NAME–VALUE PAIRS (optional)
%   'BaseDir'      : base folder (default: pwd)
%   'FTPath'       : FieldTrip root (default: '')
%   'Fs'           : sampling rate in Hz (default: 600)
%   'MinTrials'    : minimum trials required (default: 5)
%   'NW'           : DPSS time-bandwidth product (default: 3)
%   'K'            : number of tapers (default: [])  (if empty, FieldTrip picks)
%   'Pad'          : padding for FFT (default: 'nextpow2')
%   'Debias'       : debias CSD (default: false) (kept for compatibility; not always used by FT)
%
% OUTPUT FILES
%   <BaseDir>/<subj>/psi_<cond>_<freqLabel>.mat
%     - psi_mat [nROI x nROI]   (single band summary; see "PSI reduction" below)
%     - freqs   [1 x nFreq]     (frequency axis used internally)
%     - ROIs, artLog, meta
%   <BaseDir>/<subj>/psi_<cond>_<freqLabel>_summary.txt
%
% NOTES
%   - PSI is computed from the (complex) cross-spectral density across multiple
%     frequencies. Therefore we compute a spectrum over freqRange (not a single f0).
%   - We use mtmfft + DPSS to estimate CSD, then ft_connectivityanalysis(method='psi').
%   - PSI output is typically psi(spctrm) with dimord 'chan_chan_freq'. We reduce it to
%     a single matrix by averaging across frequencies within freqRange.
%   - Diagonal is set to NaN (self-coupling not meaningful).
%
% FIELDTRIP PATH HYGIENE
%   Use only ONE FieldTrip version on your MATLAB path and avoid addpath(genpath()).

% -------------------- parse args --------------------
p = inputParser;
p.addParameter('BaseDir',   pwd,                                   @(x)ischar(x)||isstring(x));
p.addParameter('FTPath',    '',                                    @(x)ischar(x)||isstring(x));
p.addParameter('Fs',        600,                                   @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinTrials', 5,                                     @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('NW',        3,                                     @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('K',         [],                                    @(x)isnumeric(x)&&isscalar(x)&&x>=1 || isempty(x));
p.addParameter('Pad',       'nextpow2',                             @(x)ischar(x)||isstring(x));
p.addParameter('Debias',    false,                                 @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

opt.BaseDir = char(opt.BaseDir);
opt.FTPath  = char(opt.FTPath);
opt.Pad     = char(opt.Pad);

% -------------------- toolboxes ---------------------
if ~isempty(opt.FTPath) && exist(opt.FTPath,'dir')
    addpath(opt.FTPath);
end
if exist('ft_defaults','file')~=2
    error('FieldTrip not found. Provide FTPath or add FieldTrip root to the MATLAB path.');
end
ft_defaults;

% -------------------- paths -------------------------
subjDir = fullfile(opt.BaseDir, subj);
infile  = fullfile(subjDir, sprintf('roi_timeseries_%s_%s.mat', cond, freqLabel));

if ~exist(infile, 'file')
    fprintf('[%s-%s-%s] ROI file not found: %s — skipping.\n', subj, cond, freqLabel, infile);
    psi_outfile = '';
    return
end

% -------------------- load --------------------------
S = load(infile); % expects roi_data, ROIs, optional artLog
if ~isfield(S,'roi_data') || isempty(S.roi_data)
    fprintf('[%s-%s-%s] roi_data missing/empty — skipping.\n', subj, cond, freqLabel);
    psi_outfile = '';
    return
end

roi_data = S.roi_data;     % [nROI x nTrials x nSamples]
ROIs     = S.ROIs;

% Optional artifact log (pass-through)
artLog = struct('nOrig', size(roi_data,2), 'nRejected', NaN, 'removed', [], 'metrics', struct());
if isfield(S,'artLog') && ~isempty(S.artLog)
    artLog = S.artLog;
end

[nROI, nTrials, nSamples] = size(roi_data);

if nTrials < opt.MinTrials
    summaryfile = fullfile(subjDir, sprintf('psi_%s_%s_summary.txt', cond, freqLabel));
    fid = fopen(summaryfile, 'w');
    if fid>0
        fprintf(fid, 'PSI for %s (%s, %s): fewer than %d valid trials — skipped.\n', ...
            subj, cond, freqLabel, opt.MinTrials);
        fclose(fid);
    end
    fprintf('[%s-%s-%s] < %d trials — skipping.\n', subj, cond, freqLabel, opt.MinTrials);
    psi_outfile = '';
    return
end

labels = {ROIs.name};
labels = labels(:);

% -------------------- build FT raw ------------------
dat = [];
dat.label    = labels;
dat.fsample  = opt.Fs;
dat.trial    = cell(1, nTrials);
dat.time     = cell(1, nTrials);

t = (0:nSamples-1) / opt.Fs;
for k = 1:nTrials
    X = squeeze(roi_data(:, k, :));  % [nROI x nSamples]
    dat.trial{k} = X;
    dat.time{k}  = t;
end
dat = ft_checkdata(dat, 'datatype','raw', 'hassampleinfo','yes');

% -------------------- frequency support --------------------
fmin = freqRange(1); fmax = freqRange(2);
if ~(isnumeric(freqRange) && numel(freqRange)==2 && fmin>0 && fmax>fmin)
    error('Invalid freqRange: must be [fmin fmax] with 0 < fmin < fmax.');
end

% We compute a frequency axis within the band.
% Keep it simple and robust: integer Hz grid within band.
foi = ceil(fmin) : 1 : floor(fmax);
if numel(foi) < 3
    % PSI needs multiple frequency points; enforce at least 3 points.
    foi = linspace(fmin, fmax, 5);
end

% -------------------- spectral estimation (CSD) --------------------
cfg_fft = [];
cfg_fft.method      = 'mtmfft';
cfg_fft.output      = 'fourier';      % keeps trials -> enables CSD computation inside FT
cfg_fft.taper       = 'dpss';
cfg_fft.foi         = foi;
cfg_fft.pad         = opt.Pad;
cfg_fft.keeptrials  = 'yes';

% DPSS settings: either specify tapsmofrq (common) or taper count.
% We use NW/K when provided to keep control explicit.
% FieldTrip can accept cfg.tapsmofrq or cfg.taper='dpss' only; NW/K are not always direct cfg fields.
% Robust approach: use tapsmofrq derived from band width, but cap to a sensible value.
bw = fmax - fmin;
cfg_fft.tapsmofrq = max(1, bw/4); % moderate smoothing; ensures stable CSD across band

freq = ft_freqanalysis(cfg_fft, dat);

% Ensure we have a cross-spectrum available for PSI:
% FieldTrip may store Fourier coefficients; connectivityanalysis will derive CSD internally.
% But PSI requires CSD between all channel pairs; ft_connectivityanalysis(method='psi') handles this.

% -------------------- PSI --------------------------
cfg_con = [];
cfg_con.method = 'psi';

% IMPORTANT:
% - Do NOT pass cfg.channelcmb unless you truly need it; passing it wrong triggers
%   "required autocombination not found" errors.
% - Keep default behaviour (all-to-all).
conn = ft_connectivityanalysis(cfg_con, freq);

% Extract PSI field robustly (varies by FT version / config)
if isfield(conn, 'psispctrm')
    PSI = conn.psispctrm;            % [nROI x nROI x nFreq]
elseif isfield(conn, 'psi')
    PSI = conn.psi;
else
    error('PSI field not found in connectivity output (expected psispctrm or psi).');
end

% Reduce to single matrix by averaging across frequency dimension (omit NaNs)
if ndims(PSI) == 3
    psi_mat = mean(PSI, 3, 'omitnan');
else
    psi_mat = PSI;
end

% Set diagonal to NaN
psi_mat(1:size(psi_mat,1)+1:end) = NaN;

% -------------------- save --------------------------
psi_outfile = fullfile(subjDir, sprintf('psi_%s_%s.mat', cond, freqLabel));

meta = struct();
meta.subj      = subj;
meta.cond      = cond;
meta.freqLabel = freqLabel;
meta.freqRange = freqRange;
meta.foi       = foi;
meta.fs        = opt.Fs;

freqs = freq.freq; %#ok<NASGU>

save(psi_outfile, 'psi_mat', 'freqs', 'ROIs', 'artLog', 'meta', '-v7.3');

% -------------------- QC summary --------------------
summaryfile = fullfile(subjDir, sprintf('psi_%s_%s_summary.txt', cond, freqLabel));
fid = fopen(summaryfile, 'w');
if fid>0
    fprintf(fid, 'PSI summary for %s (%s), band %s [%.1f–%.1f Hz]\n', subj, cond, freqLabel, fmin, fmax);
    fprintf(fid, 'Reduction: mean across %d frequency bins (foi).\n\n', numel(freqs));

    fprintf(fid, 'ROI-wise strongest directed links (by magnitude)\n');
    fprintf(fid, 'POSITIVE (source -> target):\n');
    fprintf(fid, 'ROI\t\tTarget\t\tValue\n');
    fprintf(fid, '----------------------------------------------\n');
    for r = 1:nROI
        row = psi_mat(r,:); row(r) = NaN;
        [maxval, maxidx] = max(row, [], 'omitnan');
        if ~isnan(maxval)
            fprintf(fid, '%-12s\t%-12s\t% .4f\n', ROIs(r).name, ROIs(maxidx).name, maxval);
        else
            fprintf(fid, '%-12s\t%-12s\t% .4f\n', ROIs(r).name, 'NA', NaN);
        end
    end

    fprintf(fid, '\nNEGATIVE (source -> target):\n');
    fprintf(fid, 'ROI\t\tTarget\t\tValue\n');
    fprintf(fid, '----------------------------------------------\n');
    for r = 1:nROI
        row = psi_mat(r,:); row(r) = NaN;
        [minval, minidx] = min(row, [], 'omitnan');
        if ~isnan(minval)
            fprintf(fid, '%-12s\t%-12s\t% .4f\n', ROIs(r).name, ROIs(minidx).name, minval);
        else
            fprintf(fid, '%-12s\t%-12s\t% .4f\n', ROIs(r).name, 'NA', NaN);
        end
    end

    fprintf(fid, '\n--- Artifact log (if available) ---\n');
    nOrig = artLog.nOrig;
    nRej  = artLog.nRejected;
    pct   = 100 * double(nRej) / max(double(nOrig),1);
    if isnan(pct), pct = NaN; end
    fprintf(fid, 'Total trials: %d, removed: %d (%.1f%%)\n', nOrig, nRej, pct);

    fclose(fid);
end

fprintf('[%s-%s-%s] PSI saved: %s\n', subj, cond, freqLabel, psi_outfile);
end