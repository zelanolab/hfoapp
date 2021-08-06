function HFOAppDefaultConfig()

%% Initialize UI data
cfg = [];

% Prefiex to add to the saving file name
cfg.savename_prefix = 'HFOEvents_';

% Show xlabel as clock time or time
% 'clock' | 'seconds'
cfg.clock_type = 'clock';

% number of x ticks
cfg.nb_xticks = 4;

% ytick label rotation on time windows
cfg.time_ytick_rotation = 0;

% font color for event info and status bar
cfg.msg_fc = 'b';


%% line style and line width of the rectangle
cfg.box_linecolor = 'r';
cfg.box_linestyle = '--';
cfg.box_linewidth = 0.75;

% line color of the cursor
cfg.cursor_linecolor = 'b';
% cursor color on spectrogram plot
cfg.cursor_linecolor_freq = 'k';
cfg.cursor_linewidth = 1;

cfg.cursor_full_range = 'no';

% default event color, color of the last will be used if
% there's more than one event types for an event
cfg.all_event_color = {[1 0 0],...
    [0 0 1],...
    [0 1 0],...
    [0.5, 0.7, 0.1],...
    [0.3, 0.1, 0]};

% 'location' | 'channel' (default), how the events appear in
% the event listbox
cfg.event_sort_method = 'channel';

% detection on multiple channels
cfg.mult_chnl_detect = 'no';

% initial time window length, in seconds
cfg.win_len = 1;

% step size when the <</< and >>/> buttons are pressed
%  relative to the window length, 0-1
cfg.win_major_step = 1;
cfg.win_minor_step = 0.25;

% length of extension for adjustment of manual events, in seconds
cfg.manual_adj_win = 0.05;

% height of the cursor line and event box
cfg.cursor_ypnt = [0, 1];
cfg.clust_ypnt = [0, 0, 1, 1, 0];

cfg.fontsize = 12;
cfg.color_main = 0.05 * ones( 1, 3);
cfg.color_current_channel = 0.2 * ones( 1, 3);

% Line color between channels in main window
cfg.color_main_separator = [0.8, 0, 0];
% line color between filtered time series in main window
cfg.color_minor_separator = [0.1, 0, 0];
cfg.show_main_separator = 'yes';
cfg.show_minor_separator = 'yes';
cfg.color_current_background = 0.95 * ones( 1, 3);
cfg.show_current_background = 'yes';

% 1 x 2 cell, {label, index in bpfilt.mat}
cfg.current_channel = {};

% If there's any data point within +/- 500 ms aroudn the cusrolor meet the criteria
% cluster search time window
%search_len = 2;
cfg.autoev_search_len = 2;

% cusor must be within the cluster, or within a cenrtain
% distance from the cluster
cfg.autoev_distance_thresh = 0.5;

cfg.default_autoev_search_len = 1; % 1 s
cfg.default_autoev_distance_thresh = 0.5; % 0.5s

% merge detected events
% If all z values between events are greater than 2
% and the distance between them is less thann a threshold
% is_merge = 'no';
cfg.autoev_automerge = 'no';

% in seconds, any clusters within this distance will be mreged
% merge_durthresh = 0.05;
cfg.autoev_automerge_durthresh = 0.05;
% all samples between two clusters are greater than this threshold
cfg.autoev_automerge_zthresh = 1.49999999;
cfg.autoev_use_spec = false;

%% Time-frequency plots
% keep the spectrogram of all calculated channels
cfg.tf_keepall = 'no';

% calculate spectrogram for entire time series
% 'yes' | 'no' (current window only)
% wired showing image problem always yes maybe
cfg.tf_entire = 'yes';

% pad size for current window spectrogram calculation, in seconds
cfg.tf_padsize = 1;

% calculate spectrogram and save it to data structure
% frequency range
cfg.tf_foi_range = [40, 250];

% number of linearly spaced of frequencies
cfg.tf_foi_len = 100;

% frequency tick
cfg.tf_ytick = [];

% multiple spectrogram windows 'yes'|'no'
cfg.tf_multi_wins = 'no';

cfg.tf_autoscale = 'yes';
cfg.freq_colormap = 'jet';

cfg.entire_tf_plotted = 'yes';

cfg.thresh_inclusion = 5;
cfg.thresh_onset = 3;
cfg.thresh_nbcycle = 3;
cfg.bk_whole = 'no';
cfg.bk_current = 'no';
cfg.bk_specified = 'yes';
cfg.bk_specifiy_val = 5;


% p = fileparts( mfilename( 'fullpath'));
p = fileparts( which( 'HFOApp'));
save( fullfile( p, 'default_config'), 'cfg');
clear cfg

end

