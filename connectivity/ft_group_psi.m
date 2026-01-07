function group_outfile = ft_group_psi(subjList, cond, freqLabel, varargin)
% ft_group_psi
% Group-level aggregation of subject-wise PSI matrices (single-band files).
%
% INPUTS
%   subjList   : cell array of subject IDs (e.g., {'A0102','A0174',...})
%   cond       : condition label (e.g., 'EMG2' or 'EMG3')
%   freqLabel  : band label used in filenames (e.g., 'alpha','beta','lowGamma','highGamma')
%
% NAME–VALUE PAIRS
%   'BaseDir'     : base directory containing subject folders (default: pwd)
%   'OutDir'      : output directory (default: BaseDir)
%   'MakeFigure'  : true/false, show group heatmap (default: true)
%   'ExportCSV'   : true/false, export group mean matrix as CSV (default: true)
%   'MinSubjects' : minimum number of valid subjects required (default: 5)
%   'Verbose'     : true/false (default: true)
%
% OUTPUTS
%   group_outfile : full path to saved MAT file:
%      <OutDir>/group_psi_<cond>_<freqLabel>.mat
%
% SAVED CONTENT (MAT)
%   group_mean, group_std, all_mats, roi_labels, subj_used, subj_missing, subj_invalid
%
% Notes (repo-friendly):
% - No hard-coded lab paths
% - No subject identifiers printed unless Verbose=true
% - Robust handling of missing/invalid files
% -------------------------------------------------------------------------

% -------------------- Parse args --------------------
p = inputParser;
p.addParameter('BaseDir',     pwd,   @(x)ischar(x)||isstring(x));
p.addParameter('OutDir',      '',    @(x)ischar(x)||isstring(x));
p.addParameter('MakeFigure',  true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('ExportCSV',   true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('MinSubjects', 5,     @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Verbose',     true,  @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

opt.BaseDir = char(opt.BaseDir);
if isempty(opt.OutDir)
    opt.OutDir = opt.BaseDir;
else
    opt.OutDir = char(opt.OutDir);
end
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

assert(iscell(subjList) && ~isempty(subjList), 'subjList must be a non-empty cell array.');
cond      = char(cond);
freqLabel = char(freqLabel);

% -------------------- Collect subject data --------------------
nSubj = numel(subjList);

all_mats      = [];     % will become nROI x nROI x nValid
roi_labels    = {};
subj_used     = {};
subj_missing  = {};
subj_invalid  = {};

for i = 1:nSubj
    subj = subjList{i};
    infile = fullfile(opt.BaseDir, subj, sprintf('psi_%s_%s.mat', cond, freqLabel));

    if exist(infile,'file')~=2
        subj_missing{end+1} = subj; %#ok<AGROW>
        if opt.Verbose
            warning('Missing PSI file: %s', infile);
        end
        continue;
    end

    % Load minimally
    try
        S = load(infile, 'psi_mat', 'ROIs');
    catch ME
        subj_invalid{end+1} = subj; %#ok<AGROW>
        if opt.Verbose
            warning('Could not load %s (%s).', infile, ME.message);
        end
        continue;
    end

    if ~isfield(S,'psi_mat') || isempty(S.psi_mat) || ndims(S.psi_mat) < 2
        subj_invalid{end+1} = subj; %#ok<AGROW>
        if opt.Verbose
            warning('Invalid psi_mat in %s', infile);
        end
        continue;
    end
    if ~isfield(S,'ROIs') || ~isfield(S.ROIs,'name')
        subj_invalid{end+1} = subj; %#ok<AGROW>
        if opt.Verbose
            warning('Missing ROIs/name in %s', infile);
        end
        continue;
    end

    % On first valid subject, initialize dimensions
    if isempty(all_mats)
        nROI = numel(S.ROIs);
        roi_labels = {S.ROIs.name};

        if size(S.psi_mat,1) ~= nROI || size(S.psi_mat,2) ~= nROI
            subj_invalid{end+1} = subj; %#ok<AGROW>
            if opt.Verbose
                warning('Size mismatch in %s: psi_mat is %dx%d but ROIs has %d.', ...
                    infile, size(S.psi_mat,1), size(S.psi_mat,2), nROI);
            end
            all_mats = []; roi_labels = {};
            continue;
        end

        all_mats = nan(nROI, nROI, 0, 'double');
    else
        % Enforce consistent ROI count
        if numel(S.ROIs) ~= numel(roi_labels) || size(S.psi_mat,1) ~= numel(roi_labels)
            subj_invalid{end+1} = subj; %#ok<AGROW>
            if opt.Verbose
                warning('ROI mismatch vs. reference in %s. Skipping.', infile);
            end
            continue;
        end
        % Optional: enforce consistent label order (strict)
        this_labels = {S.ROIs.name};
        if ~isequal(this_labels, roi_labels)
            subj_invalid{end+1} = subj; %#ok<AGROW>
            if opt.Verbose
                warning('ROI label/order mismatch in %s. Skipping.', infile);
            end
            continue;
        end
    end

    % Use single band/bin: accept either ROIxROI or ROIxROIx1
    M = double(S.psi_mat(:,:,1));

    % Basic sanity: keep diag as NaN (directional PSI, self-links not informative)
    M(1:size(M,1)+1:end) = NaN;

    all_mats(:,:,end+1) = M; %#ok<AGROW>
    subj_used{end+1} = subj; %#ok<AGROW>
end

nValid = size(all_mats,3);
if opt.Verbose
    fprintf('[ft_group_psi] Valid subjects: %d | missing: %d | invalid: %d\n', ...
        nValid, numel(subj_missing), numel(subj_invalid));
end

assert(nValid >= opt.MinSubjects, ...
    'Not enough valid subjects (%d) for group aggregation (MinSubjects=%d).', nValid, opt.MinSubjects);

% -------------------- Compute stats --------------------
group_mean = mean(all_mats, 3, 'omitnan');
group_std  = std(all_mats,  0, 3, 'omitnan');

% -------------------- Optional figure --------------------
if opt.MakeFigure
    fig = figure('Color','w','Position',[100 100 900 750]);
    imagesc(group_mean);
    axis square;
    set(gca,'XTick',1:numel(roi_labels),'XTickLabel',roi_labels,'XTickLabelRotation',90);
    set(gca,'YTick',1:numel(roi_labels),'YTickLabel',roi_labels);
    title(sprintf('Group mean PSI (%s, %s), n=%d', cond, freqLabel, nValid), 'Interpreter','none');
    cb = colorbar; %#ok<NASGU>
end

% -------------------- Save outputs --------------------
group_outfile = fullfile(opt.OutDir, sprintf('group_psi_%s_%s.mat', cond, freqLabel));
save(group_outfile, ...
    'group_mean','group_std','all_mats','roi_labels', ...
    'subj_used','subj_missing','subj_invalid', ...
    'cond','freqLabel','-v7.3');

if opt.Verbose
    fprintf('[ft_group_psi] Saved: %s\n', group_outfile);
end

% -------------------- CSV export --------------------
if opt.ExportCSV
    csvfile = fullfile(opt.OutDir, sprintf('group_psi_%s_%s_mean.csv', cond, freqLabel));
    writematrix(group_mean, csvfile);
    if opt.Verbose
        fprintf('[ft_group_psi] CSV exported: %s\n', csvfile);
    end
end

% -------------------- QC text summary --------------------
qcfile = fullfile(opt.OutDir, sprintf('group_psi_%s_%s_summary.txt', cond, freqLabel));
fid = fopen(qcfile, 'w');
assert(fid>0, 'Could not open QC file for writing: %s', qcfile);

fprintf(fid, 'Group PSI summary\n');
fprintf(fid, 'Condition: %s\nBand: %s\nValid subjects: %d\n\n', cond, freqLabel, nValid);

fprintf(fid, 'ROI-wise strongest directed links (by magnitude)\n');
fprintf(fid, 'Positive PSI (source -> target):\n');
fprintf(fid, 'SourceROI\tTargetROI\tValue\n');
fprintf(fid, '----------------------------------------------\n');
for r = 1:numel(roi_labels)
    row = group_mean(r,:);
    row(r) = NaN;
    [maxval, maxidx] = max(row, [], 'omitnan');
    if ~isnan(maxval)
        fprintf(fid, '%s\t%s\t% .4f\n', roi_labels{r}, roi_labels{maxidx}, maxval);
    else
        fprintf(fid, '%s\t%s\t% .4f\n', roi_labels{r}, 'NA', NaN);
    end
end

fprintf(fid, '\nNegative PSI (source -> target):\n');
fprintf(fid, 'SourceROI\tTargetROI\tValue\n');
fprintf(fid, '----------------------------------------------\n');
for r = 1:numel(roi_labels)
    row = group_mean(r,:);
    row(r) = NaN;
    [minval, minidx] = min(row, [], 'omitnan');
    if ~isnan(minval)
        fprintf(fid, '%s\t%s\t% .4f\n', roi_labels{r}, roi_labels{minidx}, minval);
    else
        fprintf(fid, '%s\t%s\t% .4f\n', roi_labels{r}, 'NA', NaN);
    end
end

fprintf(fid, '\nMatrix (group mean; NaN on diagonal)\n');
fprintf(fid, 'Rows = source, cols = target\n\n');

% Header row
fprintf(fid, '\t');
for c = 1:numel(roi_labels)
    fprintf(fid, '%s\t', roi_labels{c});
end
fprintf(fid, '\n');

for r = 1:numel(roi_labels)
    fprintf(fid, '%s\t', roi_labels{r});
    for c = 1:numel(roi_labels)
        if isnan(group_mean(r,c))
            fprintf(fid, 'NaN\t');
        else
            fprintf(fid, '%.4f\t', group_mean(r,c));
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);

if opt.Verbose
    fprintf('[ft_group_psi] QC summary saved: %s\n', qcfile);
end

end