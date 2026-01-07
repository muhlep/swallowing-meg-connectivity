function [trialData, trl] = ft_GenTrialData(dsname, varargin)
% ft_GenTrialData
% Trial definition + preprocessing for MEG data recorded on CTF systems (.ds).
% Uses FieldTrip to define trials based on an event label (e.g., EMG2/EMG3),
% then applies optional baseline correction, line-noise DFT filtering, and
% band-pass filtering.
%
% INPUT
%   dsname : path to CTF dataset (.ds directory). If empty, a folder dialog opens.
%
% NAME–VALUE PAIRS
%   'triggerName'    : eventtype for ft_definetrial (default: 'EMG2')
%   'epoch'          : [tmin tmax] in seconds, e.g. [-0.4 0.6] (default: [])
%   'baseline'       : [t1 t2] in seconds for demean baselinewindow (default: [])
%   'channel'        : channels to include (default: 'MEG')
%   'removeLinefreq' : true/false, apply DFT notch at 50/100/150 Hz (default: false)
%   'lineFreqs'      : vector of line freqs for DFT (default: [50 100 150])
%   'padding'        : padding (s) for DFT (default: 1.0)
%   'bpFreq'         : [f1 f2] band-pass (Hz) (default: [])
%   'outFile'        : optional .mat output path (default: '')
%
% OUTPUTS
%   trialData : FieldTrip data structure after ft_preprocessing
%   trl       : cfg.trl (trial definition matrix)
%
% REQUIREMENTS
%   FieldTrip must be on the MATLAB path (ft_defaults available).
%
% EXAMPLE
%   [D, trl] = ft_GenTrialData('/path/subX.ds', ...
%       'triggerName','EMG2', ...
%       'epoch',[-0.4 0.6], ...
%       'baseline',[-0.4 0], ...
%       'removeLinefreq',true, ...
%       'bpFreq',[8 13]);
%
% Notes (repo-friendly):
% - No hard-coded lab paths
% - No subject lists embedded
% - No project-specific filenames required
% -------------------------------------------------------------------------

% ----------------------------- Defaults -----------------------------
if nargin < 1, dsname = ''; end

opt.triggerName    = 'EMG2';
opt.epoch          = [];         % [tmin tmax] in seconds (tmin can be negative)
opt.baseline       = [];         % [t1 t2] in seconds
opt.channel        = 'MEG';
opt.removeLinefreq = false;
opt.lineFreqs      = [50 100 150];
opt.padding        = 1.0;
opt.bpFreq         = [];
opt.outFile        = '';

% ----------------------------- Parse varargin -------------------------
assert(mod(numel(varargin),2)==0, 'Arguments must be name–value pairs.');

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};

    switch key
        case "triggername"
            opt.triggerName = char(val);

        case "epoch"
            validateattributes(val, {'numeric'}, {'vector','numel',2}, mfilename, 'epoch');
            opt.epoch = val(:)';

        case "baseline"
            validateattributes(val, {'numeric'}, {'vector','numel',2}, mfilename, 'baseline');
            opt.baseline = val(:)';

        case "channel"
            opt.channel = val;

        case "removelinefreq"
            opt.removeLinefreq = logical(val);

        case "linefreqs"
            validateattributes(val, {'numeric'}, {'vector','positive'}, mfilename, 'lineFreqs');
            opt.lineFreqs = val(:)';

        case "padding"
            validateattributes(val, {'numeric'}, {'scalar','>',0}, mfilename, 'padding');
            opt.padding = double(val);

        case "bpfreq"
            validateattributes(val, {'numeric'}, {'vector','numel',2,'increasing','positive'}, mfilename, 'bpFreq');
            opt.bpFreq = val(:)';

        case "outfile"
            opt.outFile = char(val);

        otherwise
            error('Unknown parameter: %s', key);
    end
end

% ----------------------------- FieldTrip check ------------------------
assert(exist('ft_defaults','file')==2, 'FieldTrip not found on the MATLAB path (ft_defaults missing).');
% Do not call ft_defaults here to avoid side-effects in batch pipelines.
% (Call it once in your main entry script.)

% ----------------------------- Dataset selection ----------------------
if isempty(dsname)
    dsname = uigetdir([], 'Select a CTF dataset (.ds directory)');
    if isequal(dsname,0), error('No dataset selected.'); end
end

assert(exist(dsname,'dir')==7, 'CTF dataset not found: %s', dsname);

% ----------------------------- Trial definition -----------------------
cfg = [];
cfg.dataset            = dsname;
cfg.continuous         = 'yes';
cfg.trialdef.eventtype = opt.triggerName;

if ~isempty(opt.epoch)
    tmin = opt.epoch(1);
    tmax = opt.epoch(2);
    assert(tmin <= 0 && tmax >= 0, 'epoch must be [tmin tmax] with tmin<=0 and tmax>=0 (seconds).');
    cfg.trialdef.prestim  = abs(tmin);
    cfg.trialdef.poststim = tmax;
end

try
    cfg = ft_definetrial(cfg);
catch ME
    error('ft_definetrial failed (eventtype=%s) for %s:\n%s', opt.triggerName, dsname, ME.message);
end

trl = cfg.trl;
fprintf('[ft_GenTrialData] Trials defined: %d | eventtype=%s\n', size(trl,1), opt.triggerName);

% ----------------------------- Preprocessing --------------------------
cfg.channel = opt.channel;

if opt.removeLinefreq
    cfg.dftfilter = 'yes';
    cfg.dftfreq   = opt.lineFreqs;
    cfg.padding   = opt.padding;
end

if ~isempty(opt.baseline)
    cfg.demean         = 'yes';
    cfg.baselinewindow = opt.baseline;
end

if ~isempty(opt.bpFreq)
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = opt.bpFreq;
end

try
    trialData = ft_preprocessing(cfg);
catch ME
    error('ft_preprocessing failed for %s:\n%s', dsname, ME.message);
end

% ----------------------------- Optional save --------------------------
if ~isempty(opt.outFile)
    outFolder = fileparts(opt.outFile);
    if ~isempty(outFolder) && ~exist(outFolder,'dir'), mkdir(outFolder); end
    save(opt.outFile, 'trialData', 'trl', '-v7.3');
    fprintf('[ft_GenTrialData] Saved: %s\n', opt.outFile);
end

end