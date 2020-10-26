function [varargout] = shr_makeTargetFiles(s, b, varargin)
%
%   This function creates .tgt and .dtp files for the Sequence Horizon
%   Reaching behavioral experiment. In order to function properly, this
%   function requires both the dataframe and the plotlib toolboxes.
%   dataframe:  https://github.com/jdiedrichsen/dataframe.git
%   plotlib:    https://github.com/nejaz1/plotlib.git
%
% INPUTS
%     s:            subject number (can be vector for multiple subjects)
%     b:            block number (can be vector for multiple blocks)
%
% OPTIONAL VARARGINS
%     n_trials:     defines how many trials per block are created
%     seq_len:      sequence length (number of items / targets)
%     seq_hor:      sequence horizon (how many targets ahead are visible)
%     x_range:      spatial range of workspace x coordinates [x_min x_max]
%     y_range:      spatial range of workspace y coordinates [y_min y_max]
%     use_grid:     option of whether to use a grid of pre-determined target
%                   locations (1 = grid locations), or not (0 = random locations)
%     grid_size:    (only matters if using grid) spcifies the resolution of the grid
%                   (i.e. how many lines are splitting the workspace to create the grid)
%     min_d:        (only matters if NOT using grid) minimum distance between successive
%                   targets in percentage of largest coord range (default is 20%)
%     max_d:        (only matters if NOT using grid) maximum distance between successive
%                   targets in percentage of largest coord range (default is 21%)
%     show_seq:     option to plot the sequence or targets for each trial of each block
%                   (use carefully! e.g., by limiting the number of trials and blocks)
%
% OUTPUTS
%     G:            output struct for the whole group of subjects [subjects
%                   x blocks x trials]
%     .tgt files    (one per block per subject) specifying the order and
%                   type of trials (spreadsheet type)
%     .dtp files    (one per block per subject) specifying information to
%                   be read in by Kinarm Simulink
%
% USAGE EXAMPLES
%     [G] = shr_makeTargetFiles(99, 1, 'n_trials',1, 'use_grid',1, 'grid_size',11, 'seq_len',10, 'show_seq',1);
%     [G] = shr_makeTargetFiles(99, 1, 'n_trials',1, 'use_grid',0, 'min_d',5, 'max_d',95, 'seq_len',10, 'show_seq',1);
%     [G] = shr_makeTargetFiles([1:20], [1:12], 'use_grid',0, 'min_d',5, 'max_d',95, 'show_seq',0);
%
% --
% gariani@uwo.ca - 2020.02.24

%% make target files for all subjects and blocks per subject
G = struct();
for s = s
    S = struct();
    for b = b
        fprintf(1, '\nsubj: %d   block: %d\n', s, b);
        [B] = shr_target(s, b, varargin{:}); % B=block (all trials)
        S = addstruct(S, B); % S=subject (all blocks)
    end
    G = addstruct(G, S); % G=group (all subjects)
end
varargout{1} = G;
end

function [varargout] = shr_target(s, b, varargin)
% function [varargout] = shr_target(s, b, varargin)
% This function generates .tgt files, one per block
%
% inputs: vector of subject numbers (s), vector of block numbers (b)
% output: saved filename (fn), block structure (B)

%% define target folder
%tgtDir = '/Users/jonathanamichaels/Dropbox/KINARM/Continuous_Reach/tgt'; %save tgt files in the right folder
%if ~exist(tgtDir,'dir'); mkdir(tgtDir); end %create target folder if it doesn't already exist
dtpDir = 'C:\Users\User\Documents\Dexterit-E 3.8 Tasks\Continuous_Reach\'; %save tgt files in the right folder
if ~exist(dtpDir,'dir'); mkdir(dtpDir); end %create target folder if it doesn't already exist

%% default experimental details and varargin options
%------------------------------------------------------------------------------------------------------------------------------------
%%% IF YOU CHANGE THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADJUST first entry of <blocktable> in shr_template.dtp %%%%%%%%%%%%%%%
% !! n_trials must be a multiple of seq_hor to have a balanced number of
% horizon levels per block (e.g., if seq_hor=4, n_trials=4,8,12,16, ...)
n_trials = 2000; % how many trials in a block?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------------------------------------------------------------------
seq_length = 1; % how much lookahead?
x_range = [-4.5 4.5]; % what's the range of x coords?
y_range = [-4 4]; % what's the range of y coords?
reward_range = [150 300];
min_d = 3; % what's the minumum distance allowed between successive targets?
max_d = 8;
%%% These have to be tied to the logical radius of the targets and the minimum distance in order to avoid perturbations being assigned inside a target
min_pdist = 0.27;
max_pdist = 0.73;
%%%
load_mag = 0.3;
pert_prob = 0;
reward_prob = 1;
show_seq = 0; % do you want to plot the sequence for each trial in this block? yes (1), or no (0)
vararginoptions(varargin, {'n_trials', 'seq_length', 'x_range', 'y_range', ...
    'min_d', 'max_d', 'pert_prob', 'reward_prob', 'reward_range', 'show_seq'});

%% fill in dataframe structure B for this block
B = struct(); % initialize block structure
trial_mat = []; % initialize empty trial list
target_mat = [];
load_mat = [];
for t = 1 : n_trials
    % pick a sequence of targets with respective x,y coordinates
    % either randomly within the predefined range (some constraints apply): use_grid == 0
    % or randomly but among a predefined grid ofpotential targets: use_grid == 1
    %[x_coord, y_coord, dist, tgt_num] = pick_nrand_targs(seq_len, 'x_range',x_range, 'y_range',y_range, 'use_grid',use_grid, 'grid_size',grid_size, 'min_d',min_d, 'max_d',max_d);
    
    x = unifrnd(min(x_range), max(x_range));
    y = unifrnd(min(y_range), max(y_range));
    
    if t > 1
        last_target = target_mat(t-1,1:2);
        d = sqrt(sum((last_target - [x y]).^2));
        while d < min_d || d > max_d
            x = unifrnd(min(x_range), max(x_range));
            y = unifrnd(min(y_range), max(y_range));
            d = sqrt(sum((last_target - [x y]).^2));
        end
        if t > 2
            last_target2 = target_mat(t-2,1:2);
            d2 = sqrt(sum((last_target2 - [x y]).^2));
            while d < min_d || d > max_d || d2 < min_d || d2 > max_d
                x = unifrnd(min(x_range), max(x_range));
                y = unifrnd(min(y_range), max(y_range));
                d = sqrt(sum((last_target - [x y]).^2));
                d2 = sqrt(sum((last_target2 - [x y]).^2));
            end
        end
    end
    
    %% TARGET TABLE
    x_pos = x;
    y_pos = y;
    visual_rad = 0.8;
    logical_rad = 0.8;
    color = 9999;
    reached_color = 99999;
    
    this_target = [x_pos, y_pos, visual_rad, logical_rad, color, reached_color];
    target_mat = [target_mat; this_target];
    
    %% TP TABLE
    target_1 = t;
    target_2 = t+1;
    target_3 = t+2;
    if t == n_trials-1
        target_3 = 1;
    end
    if t == n_trials
        target_2 = 1;
        target_3 = 2;
    end
    if rand < pert_prob && t > 1
        perturb_dist = unifrnd(min_pdist, max_pdist) * d;
    else
        perturb_dist = -1;
    end
    perturb_dist = perturb_dist * 0.01; % convert to meters
    dwell_time = 200;
    if rand < reward_prob
        reward = round(unifrnd(reward_range(1), reward_range(2)));
        reward_time_l = reward;
        reward_time_u = reward;
    else
        reward_time_l = 0;
        reward_time_u = 0;
    end
    
    load_row = t;
        
    this_trial = [target_1, target_2, target_3, seq_length, perturb_dist, dwell_time, reward_time_l, reward_time_u, load_row];
    trial_mat = [trial_mat; this_trial];
    
    %% LOAD TABLE
    dir = unifrnd(0, 2*pi);
    y_load = sin(dir) * load_mag;
    x_load = cos(dir) * load_mag;
    ramp_dur = 0;
    load_dur = 100;
    this_load = [x_load, y_load, ramp_dur, load_dur];
    load_mat = [load_mat; this_load];

    if show_seq == 1
        %------------------------------------------------------------------------------------------------------------------------------------
        % sanity check: how does the sequence actually look?
   %     plot(x_coord',y_coord', 'color',[.88 .88 .88], 'linewidth',3); grid on; hold on;
   %     plt.scatter(x_coord',y_coord', 'split',tgt_num', 'regression','none', 'label',tgt_num);
   %     axis image; xlim(x_range); ylim(y_range); title(sprintf('Horizon: %d', sh(t))); hold off;
   %     set(gcf,'WindowStyle','docked');
        %------------------------------------------------------------------------------------------------------------------------------------
    end
end

%% save structure B as a target file (.tgt) and return output data structure B
outfname = sprintf('Continuous_Reach_s%02d_b%02d', s, b);
%dsave(fullfile(tgtDir, sprintf('%s.tgt', outfname)), B);

%% export and save target files in right format for exoskeleton (.dtp)
[dtp] = convert2dtp({trial_mat, target_mat, load_mat}, {'tptable', 'targettable', 'loadtable'});
xmlwrite(fullfile(dtpDir, sprintf('%s.dtp', outfname)), dtp);

%% return output structure B
B.fn = outfname;
varargout{1} = B;
end

function [dtp] = convert2dtp(allmat, tag_names)
% read in template .dtp file
dtp = xmlread('shr_template.dtp');
for i = 1:length(allmat)
    mat = allmat{i};
    mat_size = size(mat);
    tp_table = '[';
    newline = char(10);
    for row=1:mat_size(1)
        tp_table = [tp_table '['];
        for col=1:mat_size(2)
            thisData = mat(row, col);
            if thisData == 9999
                thisData = '"255255255"';
            elseif thisData == 99999
                thisData = '"000255000"';
            else
                thisData = num2str(thisData);
            end
            if col==mat_size(2)
                tp_table = [tp_table thisData];
            else
                tp_table = [tp_table thisData ', '];
            end
        end
        if row==mat_size(1)
            tp_table = [tp_table ']'];
        else
            tp_table = [tp_table '],' newline];
        end
    end
    tp_table = [tp_table ']'];
    
    % replace tp node with new tp_table info
    tp_node = dtp.getElementsByTagName(tag_names{i}).item(0).item(0);
    tp_node.set('Data', tp_table);
end
end