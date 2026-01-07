function run_all_steps(varargin)
% run_all_steps
% Minimal runner for the MEG ROI-connectivity pipeline (FieldTrip-based).
% Runs (per subject): preprocessing/trials -> ROI time series -> wPLI+PSI
% Then (group): group means and optional CBPT.
%
% IMPORTANT
% - Keep only ONE FieldTrip on the MATLAB path (no addpath(genpath(...))).
% - This runner is intentionally generic (no site-specific absolute paths).
%
% REQUIRED FUNCTIONS (same repo)
%   ft_GenTrialData.m
%   ft_CleanTrials.m
%   ft_Beamformer_extractROI.m
%   ft_wPLI_singleSub.m
%   ft_PSI_singleSub.m
%   ft_group_wPLI.m
%   ft_group_psi.m
%   ft_CBPT.m  (optional; if you keep CBPT separate, comment block below)
%
% EXAMPLE
%   run_all_steps( ...
%     'BaseDir','/path/to/Connect', ...
%     'FTPath','/path/to/fieldtrip-20231220', ...
%     'Subjects', {'A0102','A0174'}, ...
%     'DoCBPT', false);

% -------------------- user options --------------------
p = inputParser;
p.addParameter('BaseDir',  pwd, @(x)ischar(x)||isstring(x));
p.addParameter('FTPath',   '',  @(x)ischar(x)||isstring(x));
p.addParameter('SPM12Path','',  @(x)ischar(x)||isstring(x)); % only if your beamformer needs it
p.addParameter('Subjects', {},  @(x)iscell(x) || isstring(x));
p.addParameter('Conds',    {'EMG2','EMG3'}, @(x)iscell(x) || isstring(x));

% Frequency bands (labels must match your filenames & ROI time series naming)
p.addParameter('Bands', struct( ...
    'theta',    [4 8], ...
    'alpha',    [8 13], ...
    'beta',     [13 30], ...
    'lowGamma', [30 60], ...
    'highGamma',[60 80]), @(x)isstruct(x));

% Epochs (seconds) for ROI extraction (can differ per condition if you want)
p.addParameter('EpochByCond', struct( ...
    'EMG2', [-0.4 0.6], ...
    'EMG3', [0 1.0]), @(x)isstruct(x));

p.addParameter('Fs',        600, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinTrials', 5,   @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% switches
p.addParameter('DoROI',   true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('DoWPLI',  true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('DoPSI',   true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('DoGroup', true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('DoCBPT',  false, @(x)islogical(x)||ismember(x,[0 1])); % optional

% misc
p.addParameter('LineDFT', false, @(x)islogical(x)||ismember(x,[0 1])); % used in ROI extraction call
p.parse(varargin{:});
opt = p.Results;

opt.BaseDir   = char(opt.BaseDir);
opt.FTPath    = char(opt.FTPath);
opt.SPM12Path = char(opt.SPM12Path);

% -------------------- FieldTrip init (single version) --------------------
if ~isempty(opt.FTPath) && exist(opt.FTPath,'dir')
    addpath(opt.FTPath);
end
assert(exist('ft_defaults','file')==2, 'FieldTrip not found. Set FTPath correctly.');
ft_defaults;

% -------------------- subject list --------------------
subs = opt.Subjects;
if isempty(subs)
    % fallback: infer from BaseDir folders "A*"
    d = dir(fullfile(opt.BaseDir, 'A*'));
    d = d([d.isdir]);
    subs = {d.name};
end
subs = cellstr(string(subs));
fprintf('Runner: %d subjects\n', numel(subs));

conds = cellstr(string(opt.Conds));
bandNames = fieldnames(opt.Bands);

% -------------------- per-subject loop --------------------
for isub = 1:numel(subs)
    subj = subs{isub};
    fprintf('\n==================== %s (%d/%d) ====================\n', subj, isub, numel(subs));

    for ic = 1:numel(conds)
        cond = conds{ic};

        % epoch for this condition
        if isfield(opt.EpochByCond, cond)
            epoch = opt.EpochByCond.(cond);
        else
            epoch = [-1 1];
        end

        for ib = 1:numel(bandNames)
            freqLabel = bandNames{ib};
            freqBand  = opt.Bands.(freqLabel);

            fprintf('\n--- %s | %s | %s [%.1f-%.1f Hz] ---\n', subj, cond, freqLabel, freqBand(1), freqBand(2));

            % (1) ROI extraction (beamformer + ROI time series)
            if opt.DoROI
                try
                    ft_Beamformer_extractROI(subj, cond, freqBand, freqLabel, ...
                        'BaseDir',   opt.BaseDir, ...
                        'FTPath',    opt.FTPath, ...
                        'SPM12Path', opt.SPM12Path, ...
                        'Epoch',     epoch, ...
                        'MinTrials', opt.MinTrials, ...
                        'Channel',   'MEG', ...
                        'LineDFT',   opt.LineDFT);
                catch ME
                    fprintf(2, '[ROI ERROR] %s | %s | %s: %s\n', subj, cond, freqLabel, ME.message);
                end
            end

            % (2) wPLI
            if opt.DoWPLI
                try
                    ft_wPLI_singleSub(subj, cond, freqLabel, freqBand, ...
                        'BaseDir',   opt.BaseDir, ...
                        'FTPath',    opt.FTPath, ...
                        'Fs',        opt.Fs, ...
                        'MinTrials', opt.MinTrials);
                catch ME
                    fprintf(2, '[wPLI ERROR] %s | %s | %s: %s\n', subj, cond, freqLabel, ME.message);
                end
            end

            % (3) PSI
            if opt.DoPSI
                try
                    ft_PSI_singleSub(subj, cond, freqLabel, freqBand, ...
                        'BaseDir',   opt.BaseDir, ...
                        'FTPath',    opt.FTPath, ...
                        'Fs',        opt.Fs, ...
                        'MinTrials', opt.MinTrials);
                catch ME
                    fprintf(2, '[PSI ERROR] %s | %s | %s: %s\n', subj, cond, freqLabel, ME.message);
                end
            end

        end % bands
    end % conds
end % subs

% -------------------- group-level summaries --------------------
if opt.DoGroup
    fprintf('\n==================== GROUP LEVEL ====================\n');
    for ic = 1:numel(conds)
        cond = conds{ic};
        for ib = 1:numel(bandNames)
            freqLabel = bandNames{ib};

            % group mean wPLI
            if opt.DoWPLI
                try
                    ft_group_wPLI(subs, cond, freqLabel, 'BaseDir', opt.BaseDir);
                catch ME
                    fprintf(2, '[GROUP wPLI ERROR] %s | %s: %s\n', cond, freqLabel, ME.message);
                end
            end

            % group mean PSI
            if opt.DoPSI
                try
                    ft_group_psi(subs, cond, freqLabel, 'BaseDir', opt.BaseDir);
                catch ME
                    fprintf(2, '[GROUP PSI ERROR] %s | %s: %s\n', cond, freqLabel, ME.message);
                end
            end
        end
    end
end

% -------------------- optional CBPT --------------------
% If you keep CBPT as a separate “statistics” stage, this block is optional.
if opt.DoCBPT
    fprintf('\n==================== CBPT (optional) ====================\n');
    try
        % Example signature if you wrap CBPT into ft_CBPT(...)
        % ft_CBPT(subs, opt.BaseDir, 'wpli', bandNames, 'EMG2', 'EMG3');
        %
        % If your ft_CBPT takes different inputs, adapt here.
        ft_CBPT(subs, 'BaseDir', opt.BaseDir, 'Bands', bandNames);
    catch ME
        fprintf(2, '[CBPT ERROR] %s\n', ME.message);
    end
end

fprintf('\n✅ Runner finished.\n');
end