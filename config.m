
% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;

% Define quad and path to *.bbl.csv file
log_folder = 'logs';
flight_folder = '20250907';

quad = 'aosmini';
log_name = '20250907_aosmini_00.bbl.csv';

% quad = 'apex5';
% log_name = '20250907_apex5_00.bbl.csv';

% quad = 'flipmini';
% log_name = '20250908_flipmini_00.bbl.csv';

file_path = fullfile(log_folder, flight_folder, log_name);

% Evaluation parameters
do_compensate_iterm  = false;
do_show_spec_figures = true;
do_insert_legends    = true;

multp_fig_nr = ind_ax;

% Defines
set(cstprefs.tbxprefs, 'MagnitudeUnits', 'abs');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'UnwrapPhase', 'Off');
set(cstprefs.tbxprefs, 'Grid', 'On');

linewidth = 1.2;
set(0, 'defaultAxesColorOrder', get_my_colors);
pos_bode = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ... % this is a bit hacky
            0.1514, 0.1100    , 0.7536, 0.1917    ];

% Bodeoptions
opt = bodeoptions('cstprefs');