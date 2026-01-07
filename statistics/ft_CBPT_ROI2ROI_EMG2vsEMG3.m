function ft_CBPT_ROI2ROI_EMG2vsEMG3(method, varargin)
% ft_CBPT_ROI2ROI_EMG2vsEMG3
% Cluster-based permutation test (CBPT) for ROI×ROI connectivity (EMG2 vs EMG3),
% for either wPLI (undirected) or PSI (directed), across arbitrary frequency labels.
%
% This script is written to be NN-/repo-friendly:
% - no hard-coded institute paths
% - no personal identifiers in comments
% - configurable via name–value pairs
%
% EXPECTED INPUT FILES PER SUBJECT
%   <InDir>/<subj>/<method>_EMG2_<freq>.mat  contains: <method>_mat, ROIs (with .name)
%   <InDir>/<subj>/<method>_EMG3_<freq>.mat  contains: <method>_mat, ROIs (with .name)
%
% REQUIRED TOOLBOX
%   FieldTrip (tested with 2023-12-20)
%
% USAGE EXAMPLES
%   ft_CBPT_ROI2ROI_EMG2vsEMG3('wpli', 'InDir','/path/to/Connect', 'OutDir','/path/to/out');
%   ft_CBPT_ROI2ROI_EMG2vsEMG3('psi',  'SubjectList', {'S01','S02'}, 'FreqLabels', {'alpha','beta'});
%
% OUTPUTS
%   - .mat per band: stat + bookkeeping
%   - optional .xlsx cluster summary (if requested)
%
% NOTES
%   - This implements an ROI×ROI "graph space" test using ft_freqstatistics on
%     artificial 'freq-like' data (powspctrm with dimord chan_chan_freq).
%   - Neighbour structure is not used (cfg.neighbours = []).
%
% -------------------------------------------------------------------------

% ----------------------------- Parse args -----------------------------
p = inputParser;
p.FunctionName = mfilename;

validMethod = @(s) ischar(s) || isstring(s);
p.addRequired('method', validMethod);

p.addParameter('InDir',           pwd, @(x)ischar(x)||isstring(x));
p.addParameter('OutDir',          fullfile(pwd,'CBPT_results'), @(x)ischar(x)||isstring(x));

% Provide either SubjectListFile OR SubjectList
p.addParameter('SubjectListFile', '', @(x)ischar(x)||isstring(x));
p.addParameter('SubjectList',     {}, @(x)iscell(x) || isstring(x));

p.addParameter('FreqLabels',      {'theta','alpha','beta','lowGamma','highGamma'}, @(x)iscell(x)||isstring(x));

p.addParameter('NPermutations',   5000, @(x)isnumeric(x)&&isscalar(x)&&x>=100);
p.addParameter('AlphaCluster',    0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('AlphaFinal',      0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);

% Robustness options
p.addParameter('MinSubjects',     8, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('AllowNaN',        true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('MaskDiagonal',    true, @(x)islogical(x)||ismember(x,[0 1]));   % usually desired
p.addParameter('SymmetrizeWPLI',  true, @(x)islogical(x)||ismember(x,[0 1]));   % wPLI should be undirected

% Export options
p.addParameter('WriteClusterSummaryXLSX', false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('ClusterSummaryFile', '', @(x)ischar(x)||isstring(x));

% FieldTrip
p.addParameter('FTPath', '', @(x)ischar(x)||isstring(x)); % addpath(FTPath) if provided

p.parse(method, varargin{:});
opt = p.Results;

method = lower(string(opt.method));
assert(method=="wpli" || method=="psi", 'method must be ''wpli'' or ''psi''.');

opt.InDir  = char(opt.InDir);
opt.OutDir = char(opt.OutDir);
opt.FTPath = char(opt.FTPath);

% ----------------------------- FieldTrip init --------------------------
if ~isempty(opt.FTPath) && exist(opt.FTPath,'dir')==7
    addpath(opt.FTPath);
end
assert(exist('ft_defaults','file')==2, 'FieldTrip not found on MATLAB path. Add FieldTrip root or set FTPath.');
ft_defaults;

% ----------------------------- Subjects --------------------------------
subs = {};

if ~isempty(opt.SubjectList) && (iscell(opt.SubjectList) || isstring(opt.SubjectList))
    subs = cellstr(opt.SubjectList);
elseif ~isempty(opt.SubjectListFile)
    f = char(opt.SubjectListFile);
    assert(exist(f,'file')==2, 'Subject list file not found: %s', f);
    subs = strtrim(splitlines(fileread(f)));
    subs = subs(~cellfun(@isempty, subs));
else
    % fallback: list folders in InDir starting with 'A' or 'S' (common patterns)
    d = dir(opt.InDir);
    subs = {d([d.isdir] & ~startsWith({d.name},'.')).name};
end

assert(~isempty(subs), 'No subjects found. Provide SubjectList or SubjectListFile, or ensure InDir contains subject folders.');
fprintf('=== CBPT ROI×ROI EMG2 vs EMG3 | method=%s ===\n', upper(method));
fprintf('Subjects: n=%d\n', numel(subs));

% ----------------------------- Output dir ------------------------------
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

freqLabels = cellstr(opt.FreqLabels);

% Optional cluster summary output file
if opt.WriteClusterSummaryXLSX
    if ~isempty(opt.ClusterSummaryFile)
        sumfile = char(opt.ClusterSummaryFile);
    else
        sumfile = fullfile(opt.OutDir, sprintf('CBPT_ClusterSummary_%s.xlsx', upper(method)));
    end
    % init rows
    sumRows = {};
end

% ----------------------------- Main loop -------------------------------
for f = 1:numel(freqLabels)
    freqLabel = char(freqLabels{f});
    fprintf('\n--- Band: %s ---\n', freqLabel);

    group1 = {}; group2 = {}; subjIDs = {};
    nGood  = 0;

    for i = 1:numel(subs)
        subj = subs{i};

        file_EMG2 = fullfile(opt.InDir, subj, sprintf('%s_EMG2_%s.mat', method, freqLabel));
        file_EMG3 = fullfile(opt.InDir, subj, sprintf('%s_EMG3_%s.mat', method, freqLabel));

        if ~(exist(file_EMG2,'file')==2 && exist(file_EMG3,'file')==2)
            continue;
        end

        S2 = load(file_EMG2);
        S3 = load(file_EMG3);

        fieldName = sprintf('%s_mat', method);
        if ~isfield(S2, fieldName) || ~isfield(S3, fieldName)
            continue;
        end
        if ~isfield(S2,'ROIs') || ~isstruct(S2.ROIs) || ~isfield(S2.ROIs,'name')
            continue;
        end

        mat2 = S2.(fieldName);
        mat3 = S3.(fieldName);

        if ndims(mat2) > 2, mat2 = mat2(:,:,1); end
        if ndims(mat3) > 2, mat3 = mat3(:,:,1); end

        if any(size(mat2,1:2) ~= size(mat3,1:2))
            continue;
        end

        labels = {S2.ROIs.name};
        nChan  = numel(labels);
        if size(mat2,1) ~= nChan || size(mat2,2) ~= nChan
            continue;
        end

        % Clean/sanitize matrices
        mat2 = double(mat2);
        mat3 = double(mat3);

        if opt.MaskDiagonal
            mat2(1:nChan+1:end) = 0;
            mat3(1:nChan+1:end) = 0;
        end

        if method=="wpli"
            % wPLI should be non-negative undirected
            mat2 = abs(mat2); mat3 = abs(mat3);
            if opt.SymmetrizeWPLI
                mat2 = (mat2 + mat2.')/2;
                mat3 = (mat3 + mat3.')/2;
            end
        end

        if ~opt.AllowNaN
            if any(~isfinite(mat2(:))) || any(~isfinite(mat3(:)))
                continue;
            end
        end

        % Build FieldTrip "freq-like" structs for ROI×ROI CBPT:
        % powspctrm [chan x chan x freq] dimord='chan_chan_freq'
        dummy2           = [];
        dummy2.label     = labels;
        dummy2.dimord    = 'chan_chan_freq';
        dummy2.freq      = 1; % single-foi placeholder
        dummy2.powspctrm = reshape(mat2, nChan, nChan, 1);

        dummy3           = dummy2;
        dummy3.powspctrm = reshape(mat3, nChan, nChan, 1);

        group1{end+1}  = dummy2; %#ok<AGROW>
        group2{end+1}  = dummy3; %#ok<AGROW>
        subjIDs{end+1} = subj;   %#ok<AGROW>
        nGood = nGood + 1;
    end

    fprintf('Usable subjects for %s/%s: n=%d\n', upper(method), freqLabel, nGood);

    if nGood < opt.MinSubjects
        warning('Skipping %s/%s: only %d usable subjects (MinSubjects=%d).', upper(method), freqLabel, nGood, opt.MinSubjects);
        continue;
    end

    % ------------------------- CBPT config -----------------------------
    cfg = [];
    cfg.channel          = 'all';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';

    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = opt.AlphaCluster;
    cfg.clusterstatistic = 'maxsize';   % robust/simple; alternatives: 'maxsum'
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = opt.AlphaFinal;
    cfg.numrandomization = opt.NPermutations;

    % no neighbourhood in ROI×ROI space
    cfg.neighbours       = [];

    % Paired design: subject as uvar, condition as ivar
    n = numel(group1);
    cfg.design = [1:n  1:n;  ones(1,n)  2*ones(1,n)];
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    fprintf('Running CBPT (n=%d, permutations=%d)...\n', n, opt.NPermutations);
    stat = ft_freqstatistics(cfg, group1{:}, group2{:});

    % ------------------------- Save .mat -------------------------------
    outfile = fullfile(opt.OutDir, sprintf('CBPT_%s_%s_EMG2vsEMG3.mat', upper(method), freqLabel));
    save(outfile, 'stat', 'subjIDs', 'method', 'freqLabel', 'cfg', '-v7.3');
    fprintf('Saved: %s\n', outfile);

    % ------------------------- Console summary -------------------------
    report_clusters_to_console(stat);

    % ------------------------- Optional XLSX summary -------------------
    if opt.WriteClusterSummaryXLSX
        % Extract clusters (pos/neg) and list ROI pairs that belong to clusters
        rowsBand = extract_cluster_summary_rows(stat, freqLabel, method);
        sumRows  = [sumRows; rowsBand]; %#ok<AGROW>
    end
end

% ----------------------------- Write summary ---------------------------
if opt.WriteClusterSummaryXLSX
    if isempty(sumRows)
        warning('No cluster summary rows collected. Nothing to write.');
    else
        varNames = {'Method','Frequency','Contrast','ClusterType','ClusterIdx','ClusterP', ...
                    'ROI1','ROI2'};
        Tsum = cell2table(sumRows, 'VariableNames', varNames);
        writetable(Tsum, sumfile, 'FileType','spreadsheet');
        fprintf('\n✅ Cluster summary written:\n  %s\n', sumfile);
    end
end

end % function

% ======================================================================
% Local helper functions
% ======================================================================

function report_clusters_to_console(stat)
    % Positive clusters
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        pvals = arrayfun(@(c) c.prob, stat.posclusters);
        fprintf('  Positive clusters: %d | min p = %.4f\n', numel(pvals), min(pvals));
    else
        fprintf('  Positive clusters: none\n');
    end
    % Negative clusters
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        pvals = arrayfun(@(c) c.prob, stat.negclusters);
        fprintf('  Negative clusters: %d | min p = %.4f\n', numel(pvals), min(pvals));
    else
        fprintf('  Negative clusters: none\n');
    end
end

function rows = extract_cluster_summary_rows(stat, freqLabel, method)
% Returns rows:
% {Method, Frequency, Contrast, ClusterType, ClusterIdx, ClusterP, ROI1, ROI2}
%
% Uses stat.posclusterslabelmat / stat.negclusterslabelmat and stat.label.

    rows = {};
    if ~isfield(stat,'label') || isempty(stat.label)
        return;
    end
    labels = stat.label;
    nChan  = numel(labels);

    % Helper to convert cluster label matrix to ROI pairs
    function addRowsFor(typeStr, clustStruct, clustLabelMat)
        if isempty(clustStruct) || ~isfield(stat, clustLabelMat), return; end
        L = stat.(clustLabelMat);  % [chan x chan x freq] or [chan x chan]
        if ndims(L)>2, L = L(:,:,1); end

        for ci = 1:numel(clustStruct)
            p = clustStruct(ci).prob;
            if ~isfinite(p), continue; end

            mask = (L == ci);
            if ~any(mask(:)), continue; end

            [ii,jj] = find(mask);
            for k = 1:numel(ii)
                r1 = labels{ii(k)};
                r2 = labels{jj(k)};
                if ii(k)==jj(k), continue; end
                rows(end+1,:) = {upper(string(method)), string(freqLabel), "EMG2 vs EMG3", ...
                                 string(typeStr), ci, p, string(r1), string(r2)}; %#ok<AGROW>
            end
        end
    end

    if isfield(stat,'posclusters') && ~isempty(stat.posclusters) && isfield(stat,'posclusterslabelmat')
        addRowsFor("pos", stat.posclusters, 'posclusterslabelmat');
    end
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters) && isfield(stat,'negclusterslabelmat')
        addRowsFor("neg", stat.negclusters, 'negclusterslabelmat');
    end
end