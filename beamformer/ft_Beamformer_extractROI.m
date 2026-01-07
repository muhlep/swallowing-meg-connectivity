function ft_Beamformer_extractROI(subj, cond, freqBand, freq_label, varargin)
% ft_Beamformer_extractROI
% LCMV beamforming on preprocessed trials and extraction of ROI time series
% from a warped MNI sourcemodel (10 mm). Includes a QC gate: requires >=MinTrials
% after cleaning; otherwise writes an empty .mat and a QC report.
%
% INPUTS
%   subj        : subject ID string (used to locate subject folder)
%   cond        : trigger/condition label (e.g., 'EMG2' or 'EMG3')
%   freqBand    : [f_low f_high] band-pass (Hz) for preprocessing
%   freq_label  : short label for files (e.g., 'alpha','beta','lowGamma')
%
% NAME–VALUE PAIRS (optional)
%   'BaseDir'    : base directory that contains subject folders (default: pwd)
%   'Dataset'    : full path to the CTF .ds dataset (overrides default pattern)
%   'FTPath'     : FieldTrip root folder (default: '') -> ft_defaults must be reachable
%   'SPM12Path'  : SPM12 root folder (default: '')
%   'Epoch'      : [tmin tmax] in seconds for trial extraction (default: [-1 1])
%   'MinTrials'  : minimum trials required post-cleaning (default: 5)
%   'Channel'    : channel selection (default: 'MEG')
%   'LineDFT'    : logical, remove 50/100/150 Hz via DFT (default: false)
%   'Headmodel'  : path to headmodel file (default: <BaseDir>/<subj>/.../singleShell.mat)
%   'Grid'       : path to warped sourcemodel file (default: <BaseDir>/<subj>/...10mm.mat)
%   'ROIs'       : path to ROI file (default: <BaseDir>/<subj>/ROIs_onWarpedGrid.mat)
%   'OutFile'    : override output ROI timeseries .mat path (default: subject folder)
%   'QCFile'     : override QC report .txt path (default: subject folder)
%
% OUTPUTS
%   Writes:
%     - roi_timeseries_<cond>_<freq_label>.mat
%       (roi_data [nROI x nTrials x nSamples], ROIs, roi_mom_count, artLog)
%     - qc_trials_<cond>_<freq_label>.txt
%
% DEPENDENCIES
%   FieldTrip (ft_defaults), ft_GenTrialData.m, ft_CleanTrials.m.
%   Requires precomputed subject-specific headmodel + warped grid + ROI definition.
%
% NOTE
%   Ensure ONLY ONE FieldTrip version is on your MATLAB path.

% ----------------------------- Parse args -----------------------------
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('BaseDir',    pwd, @(x)ischar(x)||isstring(x));
p.addParameter('Dataset',    '',  @(x)ischar(x)||isstring(x));
p.addParameter('FTPath',     '',  @(x)ischar(x)||isstring(x));
p.addParameter('SPM12Path',  '',  @(x)ischar(x)||isstring(x));
p.addParameter('Epoch',      [-1 1], @(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x)));
p.addParameter('MinTrials',  5, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Channel',    'MEG', @(x)ischar(x)||isstring(x));
p.addParameter('LineDFT',    false, @(x)islogical(x)||ismember(x,[0 1]));

% Optional explicit input file overrides (portable repo)
p.addParameter('Headmodel',  '', @(x)ischar(x)||isstring(x));
p.addParameter('Grid',       '', @(x)ischar(x)||isstring(x));
p.addParameter('ROIs',       '', @(x)ischar(x)||isstring(x));

% Optional explicit output overrides
p.addParameter('OutFile',    '', @(x)ischar(x)||isstring(x));
p.addParameter('QCFile',     '', @(x)ischar(x)||isstring(x));

p.parse(varargin{:});
opt = p.Results;

% Normalize to char for fullfile compatibility
opt.BaseDir   = char(opt.BaseDir);
opt.Dataset   = char(opt.Dataset);
opt.FTPath    = char(opt.FTPath);
opt.SPM12Path = char(opt.SPM12Path);
opt.Channel   = char(opt.Channel);
opt.Headmodel = char(opt.Headmodel);
opt.Grid      = char(opt.Grid);
opt.ROIs      = char(opt.ROIs);
opt.OutFile   = char(opt.OutFile);
opt.QCFile    = char(opt.QCFile);

% ----------------------------- Toolboxes ------------------------------
if ~isempty(opt.FTPath) && exist(opt.FTPath,'dir')==7
    addpath(opt.FTPath);
end
assert(exist('ft_defaults','file')==2, ...
    'FieldTrip not found on MATLAB path. Add FieldTrip root (FTPath).');
ft_defaults;

if ~isempty(opt.SPM12Path) && exist(opt.SPM12Path,'dir')==7
    addpath(opt.SPM12Path);
end

% ----------------------------- Paths ----------------------------------
subjDir = fullfile(opt.BaseDir, subj);

% Dataset: override or fallback pattern
if ~isempty(opt.Dataset)
    dsname = opt.Dataset;
else
    % Fallback: common pattern (adjust in your own pipeline if needed)
    dsname = fullfile(subjDir, [subj '.ds']);
end
assert(exist(dsname,'dir')==7, 'CTF dataset not found: %s', dsname);

% Inputs: headmodel/grid/ROIs (override or fallback)
if ~isempty(opt.Headmodel)
    headmodel_file = opt.Headmodel;
else
    headmodel_file = fullfile(subjDir, 'volumeConductorModels', 'singleShell.mat');
end
if ~isempty(opt.Grid)
    grid_file = opt.Grid;
else
    grid_file = fullfile(subjDir, 'warped_MNI_sourcemodel3d10mm.mat');
end
if ~isempty(opt.ROIs)
    roi_file = opt.ROIs;
else
    roi_file = fullfile(subjDir, 'ROIs_onWarpedGrid.mat');
end

assert(exist(headmodel_file,'file')==2, 'Headmodel not found: %s', headmodel_file);
assert(exist(grid_file,'file')==2,     'Sourcemodel not found: %s', grid_file);
assert(exist(roi_file,'file')==2,      'ROI file not found: %s', roi_file);

% Outputs: default into subject folder unless overridden
if ~isempty(opt.QCFile)
    qcfile = opt.QCFile;
else
    qcfile = fullfile(subjDir, sprintf('qc_trials_%s_%s.txt', cond, freq_label));
end
if ~isempty(opt.OutFile)
    outfile = opt.OutFile;
else
    outfile = fullfile(subjDir, sprintf('roi_timeseries_%s_%s.mat', cond, freq_label));
end

% ----------------------------- Load headmodel -------------------------
S = load(headmodel_file);
if isfield(S,'vol')
    vol = S.vol;
elseif isfield(S,'headmodel')
    vol = S.headmodel;
else
    error('No ''vol'' or ''headmodel'' found in %s', headmodel_file);
end

% ----------------------------- Load grid & ROIs -----------------------
G = load(grid_file);
if isfield(G,'warped_sourcemodel')
    warped_sourcemodel = G.warped_sourcemodel;
elseif isfield(G,'sourcemodel')
    warped_sourcemodel = G.sourcemodel;
else
    error('Could not find warped sourcemodel variable in %s (expected ''warped_sourcemodel'' or ''sourcemodel'')', grid_file);
end

R = load(roi_file);
if isfield(R,'ROIs')
    ROIs = R.ROIs;
else
    error('ROI file %s does not contain variable ''ROIs''.', roi_file);
end

% ----------------------------- [1] Trials & cleaning -------------------
genArgs = { ...
    'triggerName', cond, ...
    'epoch',       opt.Epoch, ...
    'bpFreq',      freqBand, ...
    'channel',     opt.Channel};

if opt.LineDFT
    % robust: explicit boolean
    genArgs = [genArgs, {'removeLinefreq', true}]; %#ok<AGROW>
end

[data_cond, ~] = ft_GenTrialData(dsname, genArgs{:});

% Artifact cleaning
[cleanedData_cond, artLog] = ft_CleanTrials(data_cond);

% ----------------------------- QC gate ---------------------------------
nBefore   = numel(data_cond.trial);
nAfter    = numel(cleanedData_cond.trial);
minTrials = opt.MinTrials;

if nAfter < minTrials
    fid = fopen(qcfile,'w');
    if fid>0
        fprintf(fid, 'Subject: %s | Condition: %s | Freq: %s\n', subj, cond, freq_label);
        fprintf(fid, 'Trials before cleaning: %d\n', nBefore);
        fprintf(fid, 'Trials after cleaning : %d\n', nAfter);
        fprintf(fid, 'Abort: fewer than %d trials after artifact rejection.\n', minTrials);
        if isfield(artLog,'removed') && ~isempty(artLog.removed)
            fprintf(fid, 'Removed trial indices: %s\n', mat2str(artLog.removed));
        end
        fclose(fid);
    end

    roi_data      = [];
    roi_mom_count = [];
    save(outfile, 'roi_data', 'ROIs', 'roi_mom_count', 'artLog', '-v7.3');
    fprintf('[QC] %s\n', qcfile);
    fprintf('[INFO] Empty ROI file written (insufficient trials): %s\n', outfile);
    return
end

% ----------------------------- [2] Timelock ----------------------------
cfg_tl = [];
cfg_tl.covariance       = 'yes';
cfg_tl.covariancewindow = 'all';
cfg_tl.keeptrials       = 'yes';
timelock_cond = ft_timelockanalysis(cfg_tl, cleanedData_cond);

% ----------------------------- [3] Leadfield ---------------------------
cfg_lf = [];
cfg_lf.sourcemodel = warped_sourcemodel;
cfg_lf.headmodel   = vol;
cfg_lf.channel     = cleanedData_cond.label;
if isfield(cleanedData_cond,'grad'), cfg_lf.grad = cleanedData_cond.grad; end
cfg_lf.reducerank  = 2;      % typical for CTF systems
cfg_lf.normalize   = 'yes';
cfg_lf.unit        = 'mm';
leadfield = ft_prepare_leadfield(cfg_lf, cleanedData_cond);

% ----------------------------- [4] LCMV --------------------------------
cfg_bf = [];
cfg_bf.method          = 'lcmv';
cfg_bf.sourcemodel     = leadfield;
cfg_bf.headmodel       = vol;
cfg_bf.lcmv.fixedori   = 'yes';
cfg_bf.lcmv.keepmom    = 'yes';
cfg_bf.lcmv.keepfilter = 'no';
cfg_bf.keeptrials      = 'yes';
cfg_bf.rawtrial        = 'yes';
source_cond = ft_sourceanalysis(cfg_bf, timelock_cond);

% ----------------------------- [5] ROI extraction ----------------------
nROI    = numel(ROIs);
nTrials = numel(source_cond.trial);

% infer nSamples robustly from first trial with non-empty mom
example_idx = [];
if nTrials >= 1 && isfield(source_cond.trial(1),'mom') && ~isempty(source_cond.trial(1).mom)
    momCell = source_cond.trial(1).mom;
    example_idx = find(~cellfun(@isempty, momCell), 1, 'first');
end
assert(~isempty(example_idx), 'No non-empty source moments found in first trial.');

nSamples = size(source_cond.trial(1).mom{example_idx}, 2);
roi_data      = nan(nROI, nTrials, nSamples, 'double');
roi_mom_count = zeros(nROI,1);

fprintf('\n[ROI] Extracting ROI time series from warped grid (%d ROIs)...\n', nROI);

for r = 1:nROI
    if ~isfield(ROIs(r),'voxidx') || isempty(ROIs(r).voxidx)
        roi_mom_count(r) = 0;
        continue
    end

    vIdx = ROIs(r).voxidx(:);
    roi_mom_count(r) = numel(vIdx);

    for t = 1:nTrials
        moms_t = source_cond.trial(t).mom(vIdx);
        moms_t = moms_t(~cellfun(@isempty, moms_t));
        if isempty(moms_t), continue; end

        moms_mat = cellfun(@double, moms_t, 'UniformOutput', false);
        moms_mat = cat(1, moms_mat{:});              % [nVox x nSamples]
        roi_data(r, t, :) = mean(moms_mat, 1, 'omitnan');
    end
end

% ----------------------------- [6] QC report ---------------------------
fid = fopen(qcfile,'w');
if fid>0
    fprintf(fid, 'Subject: %s | Condition: %s | Freq: %s\n', subj, cond, freq_label);
    fprintf(fid, 'Trials before cleaning: %d\n', nBefore);
    fprintf(fid, 'Trials after cleaning : %d\n', nAfter);
    fprintf(fid, 'Removed trials: %d (%.1f%%)\n', nBefore-nAfter, 100*max(nBefore-nAfter,0)/max(nBefore,1));
    if isfield(artLog,'removed') && ~isempty(artLog.removed)
        fprintf(fid, 'Removed trial indices: %s\n', mat2str(artLog.removed));
    end
    fclose(fid);
end

% ----------------------------- [7] Save -------------------------------
save(outfile, 'roi_data', 'ROIs', 'roi_mom_count', 'artLog', '-v7.3');
fprintf('[SAVE] ROI time series: %s\n', outfile);
fprintf('[QC]   Trial report   : %s\n', qcfile);
fprintf('[INFO] Ready for wPLI/PSI estimation per subject.\n');
end