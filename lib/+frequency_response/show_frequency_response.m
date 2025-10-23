
function show_frequency_response(Cpi, Cd, Cpi_ana, Cd_ana, omega_bode, opt, do_insert_legends, linewidth, expand_multiple_figure_nr, multp_fig_nr)
%SHOW_FREQUENCY_RESPONSE Plots the frequency response of PI and D controllers
%
% Inputs:
%   Cpi, Cd            - measured controllers (Bode data)
%   Cpi_ana, Cd_ana    - analytical controllers (Bode data)
%   omega_bode         - frequency vector
%   opt                - structure with options (YLim, MagScale)
%   do_insert_legends  - boolean to control legend display
%   linewidth          - line width for plot
%   expand_multiple_figure_nr - function handle for figure scaling
%   multp_fig_nr       - multiplier for figure number

    % reate new figure with scaled figure number
    figure(expand_multiple_figure_nr(5, multp_fig_nr))
    
    % Bode-Plot
    bode(Cpi, Cd, Cpi_ana, Cd_ana, omega_bode, opt)
    title('Cpi, Cd')
    set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

    % add legends if requested
    if do_insert_legends
        legend('PI gemessen', 'D gemessen', 'PI analytisch', 'D analytisch')
    end
end

