function wpli_outfile = ft_wPLI_singleSub(subj, cond, freqLabel, freqRange, varargin)
% ft_wPLI_singleSub
% Compute subject-level debiased wPLI for one subject/condition/frequency band.
%
% Expects ROI time series saved by the beamformer/ROI step:
%   roi_timeseries_<cond>_<freqLabel>.mat  containing
%       roi_data : [nROI x nTrials x nSamples] (double)
%       ROIs     : struct array with field .name
%       artLog   : optional (trial cleaning log)
%
% INPUTS
%   subj       : subject ID string (used to build paths)
%   cond       : condition label (e.g., 'EMG2', 'EMG3')
%   freqLabel  : band label used in filenames (e.g., 'alpha','beta','lowGamma','highGamma')
%   freqRange  : [fmin fmax] in Hz (used for DPSS smoothing)
%
% NAME–VALUE PAIRS
%   'BaseDir'     : base folder containing subject folders (default: pwd)
%   'FTPath'      : FieldTrip root (default: '' = assume already on path)
%   'Fs'          : sampling rate in Hz (default: [] -> infer from saved time vector if present, else 600)
%   'MinTrials'   : minimum trials required (default: 5)
%   'OutFile'     : override output .mat path (default: '')
%   'Verbose'     : true/false (default: true)
%
% OUTPUT
%   wpli_outfile : path to saved <subj>/wpli_<cond>_<freqLabel>.mat
%                  (empty string if skipped)
%
% SAVED (MAT)
%   wpli_mat [nROI x nROI] (single band/bin)
%   freqs    scalar Hz (center frequency)
%   ROIs, artLog
%   meta     struct with analysis settings and provenance
%
% Notes (repo-friendly):
% - No hard-coded lab paths
% - Defensive FieldTrip output handling across versions
% - Writes a small text summary next to the MAT file
% -------------------------------------------------------------------------

% -------------------- Parse args --------------------
p = inputParser;
p.addParameter('BaseDir',    pwd,  @(x)ischar(x)||isstring(x));
p.addParameter('FTPath',     '',   @(x)ischar(x)||isstring(x));
p.addParameter('Fs',         [],   @(x)isnumeric(x)&&isscalar(x) || isempty(x));
p.addParameter('MinTrials',  5,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('OutFile',    '',   @(x)ischar(x)||isstring(x));
p.addParameter('Verbose',    true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

opt.BaseDir = char(opt.BaseDir);
opt.FTPath  = char(opt.FTPath);
if ~isempty(opt.OutFile), opt.OutFile = char(opt.OutFile); end

subj      = char(subj);
cond      = char(cond);
freqLabel = char(freqLabel);

validateattributes(freqRange, {'numeric'}, {'vector','numel',2,'finite'}, mfilename, 'freqRange');
fmin = freqRange(1); fmax = freqRange(2);
assert(fmax > fmin && fmin > 0, 'Invalid freqRange. Expect [fmin fmax] with 0 < fmin < fmax.');
f0 = mean(freqRange);
bw = fmax - fmin;

% -------------------- FieldTrip init --------------------
if ~isempty(opt.FTPath)
    if exist(opt.FTPath,'dir')~=7
        error('FTPath does not exist: %s', opt.FTPath);
    end
    addpath(opt.FTPath);
end
if exist('ft_defaults','file')~=2
    error('FieldTrip not found. Add it to your MATLAB path or pass ''FTPath''.');
end
ft_defaults;

% -------------------- Paths -------------------------
subjDir = fullfile(opt.BaseDir, subj);
infile  = fullfile(subjDir, sprintf('roi_timeseries_%s_%s.mat', cond, freqLabel));
if exist(infile,'file')~=2
    if opt.Verbose
        fprintf('[ft_wPLI_singleSub] Missing ROI input: %s (skip)\n', infile);
    end
    wpli_outfile = '';
    return;
end

if isempty(opt.OutFile)
    wpli_outfile = fullfile(subjDir, sprintf('wpli_%s_%s.mat', cond, freqLabel));
else
    wpli_outfile = opt.OutFile;
end
summaryfile = fullfile(fileparts(wpli_outfile), sprintf('wpli_%s_%s_summary.txt', cond, freqLabel));

% -------------------- Load ROI time series --------------------
S = load(infile); % expects roi_data, ROIs, optional artLog, optional roi_time/fsample

if ~isfield(S,'roi_data') || isempty(S.roi_data)
    if opt.Verbose
        fprintf('[ft_wPLI_singleSub] roi_data missing/empty in %s (skip)\n', infile);
    end
    wpli_outfile = '';
    return;
end
if ~isfield(S,'ROIs') || ~isfield(S.ROIs,'name')
    error('ROIs with .name missing in %s', infile);
end

roi_data = S.roi_data; % [nROI x nTrials x nSamples]
ROIs     = S.ROIs;

[nROI, nTrials, nSamples] = size(roi_data);
if nTrials < opt.MinTrials
    fid = fopen(summaryfile,'w');
    if fid>0
        fprintf(fid, 'wPLI skipped: fewer than %d trials.\n', opt.MinTrials);
        fprintf(fid, 'Input file: %s\n', infile);
        fprintf(fid, 'Trials available: %d\n', nTrials);
        fclose(fid);
    end
    if opt.Verbose
        fprintf('[ft_wPLI_singleSub] %s %s %s: only %d trials (<%d). Skip.\n', ...
            subj, cond, freqLabel, nTrials, opt.MinTrials);
    end
    wpli_outfile = '';
    return;
end

% artifact log (optional)
if isfield(S,'artLog') && isstruct(S.artLog)
    artLog = S.artLog;
else
    artLog = struct();
    artLog.nOrig = nTrials;
    artLog.nRejected = NaN;
    artLog.removed = [];
end

% sampling rate: prefer explicit info in file, else user value, else 600
Fs = [];
if isfield(S,'fsample') && isnumeric(S.fsample) && isscalar(S.fsample) && S.fsample>0
    Fs = S.fsample;
elseif isfield(S,'Fs') && isnumeric(S.Fs) && isscalar(S.Fs) && S.Fs>0
    Fs = S.Fs;
elseif ~isempty(opt.Fs)
    Fs = opt.Fs;
else
    Fs = 600; % fallback default used in your pipeline
end

labels = {ROIs.name};
labels = labels(:);

% -------------------- Build FT raw structure --------------------
dat = [];
dat.label   = labels;
dat.fsample = Fs;
dat.trial   = cell(1, nTrials);
dat.time    = cell(1, nTrials);

t = (0:nSamples-1) ./ Fs;
for k = 1:nTrials
    X = squeeze(roi_data(:,k,:));  % [nROI x nSamples]
    X = double(X);

    % If a trial has any all-NaN channels, keep but FT will propagate NaNs; better to guard:
    % Replace all-NaN channels with zeros to avoid complete failure, but mark via meta later.
    allNaN = all(isnan(X), 2);
    if any(allNaN)
        X(allNaN,:) = 0;
    end

    dat.trial{k} = X;
    dat.time{k}  = t;
end

dat = ft_checkdata(dat, 'datatype','raw', 'hassampleinfo','yes');

% -------------------- Frequency analysis (Fourier) --------------------
cfg_fft = [];
cfg_fft.method     = 'mtmfft';
cfg_fft.output     = 'fourier';
cfg_fft.taper      = 'dpss';
cfg_fft.foi        = f0;          % single bin at center frequency
cfg_fft.tapsmofrq  = bw/2;        % smoothing half-bandwidth
cfg_fft.pad        = 'nextpow2';
cfg_fft.keeptrials = 'yes';

freq = ft_freqanalysis(cfg_fft, dat);

% -------------------- Connectivity: debiased wPLI --------------------
cfg_con = [];
cfg_con.method  = 'wpli_debiased';
cfg_con.complex = 'absimag';

conn = ft_connectivityanalysis(cfg_con, freq);

% Robust field extraction across FT versions
C = [];
if isfield(conn,'wpli_debiasedspctrm')
    C = conn.wpli_debiasedspctrm;
elseif isfield(conn,'wplispctrm')
    C = conn.wplispctrm;
elseif isfield(conn,'wpli')
    C = conn.wpli;
elseif isfield(conn,'dwpli')
    C = conn.dwpli;
else
    % Helpful debug info
    fn = fieldnames(conn);
    error('No wPLI field found in conn. Fields were: %s', strjoin(fn(:).', ', '));
end

% Squeeze single frequency bin if present
C = squeeze(C);
if ndims(C)==3 && size(C,3)==1
    C = C(:,:,1);
end
if ~ismatrix(C) || size(C,1)~=nROI || size(C,2)~=nROI
    error('Unexpected wPLI matrix size after connectivity: got %s, expected %dx%d.', ...
        mat2str(size(C)), nROI, nROI);
end

% Make sure diagonal is NaN (self-coupling)
C(1:nROI+1:end) = NaN;

wpli_mat = double(C);
freqs    = f0;

meta = struct();
meta.subj      = subj;
meta.cond      = cond;
meta.freqLabel = freqLabel;
meta.freqRange = freqRange;
meta.foi       = f0;
meta.tapsmofrq = bw/2;
meta.fsample   = Fs;
meta.nROI      = nROI;
meta.nTrials   = nTrials;
meta.nSamples  = nSamples;
meta.input_roi_file = infile;

% -------------------- Save --------------------
save(wpli_outfile, 'wpli_mat','freqs','ROIs','artLog','meta','-v7.3');

% -------------------- QC summary text --------------------
fid = fopen(summaryfile,'w');
if fid>0
    fprintf(fid, 'Debiased wPLI summary\n');
    fprintf(fid, 'Subject:   %s\n', subj);
    fprintf(fid, 'Condition: %s\n', cond);
    fprintf(fid, 'Band:      %s  [%.1f–%.1f Hz]\n', freqLabel, fmin, fmax);
    fprintf(fid, 'Trials:    %d\n', nTrials);
    fprintf(fid, 'ROIs:      %d\n\n', nROI);

    fprintf(fid, 'ROI\tMax wPLI\tTarget ROI\n');
    fprintf(fid, '--------------------------------------\n');
    for r = 1:nROI
        row = wpli_mat(r,:);
        row(r) = NaN;
        [maxval, maxidx] = max(row, [], 'omitnan');
        if ~isnan(maxval)
            fprintf(fid, '%s\t%.4f\t%s\n', labels{r}, maxval, labels{maxidx});
        else
            fprintf(fid, '%s\tNaN\tNA\n', labels{r});
        end
    end

    fprintf(fid, '\nArtifact log (if available)\n');
    if isfield(artLog,'nOrig'),     fprintf(fid, 'nOrig: %s\n', num2str(artLog.nOrig)); end
    if isfield(artLog,'nRejected'), fprintf(fid, 'nRejected: %s\n', num2str(artLog.nRejected)); end
    if isfield(artLog,'removed') && ~isempty(artLog.removed)
        fprintf(fid, 'removed idx: %s\n', mat2str(artLog.removed));
    end

    fclose(fid);
end

if opt.Verbose
    fprintf('[ft_wPLI_singleSub] Saved: %s\n', wpli_outfile);
end

end