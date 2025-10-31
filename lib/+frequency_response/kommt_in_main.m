
% New approach with functions
%addpath('lib/+frequency_response');

results = estimate_frequency_response(data, ind, ind_ax, ind_eval, Ts_log, para, throttle_avg, Ts_cntr, data_sinarg);

% Access results via results.T, results.P, etc.