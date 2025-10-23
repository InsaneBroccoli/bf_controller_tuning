
function show_frequency_response(Cpi, Cd, Cpi_ana, Cd_ana, omega_bode, opt, do_insert_legends, linewidth, expand_multiple_figure_nr, multp_fig_nr)
%SHOW_FREQUENCY_RESPONSE Plots the frequency response of PI and D controllers
%
% Inputs:
%   Cpi, Cd            - gemessene Controller (Bode-Daten)
%   Cpi_ana, Cd_ana    - analytische Controller (Bode-Daten)
%   omega_bode         - Frequenzvektor
%   opt                - Struktur mit Optionen (YLim, MagScale)
%   do_insert_legends  - Boolean für Legende
%   linewidth          - Linienbreite
%   expand_multiple_figure_nr - Funktion handle zur Figureskalierung
%   multp_fig_nr       - Multiplikator für Figures

    % Figure erstellen
    figure(expand_multiple_figure_nr(5, multp_fig_nr))
    
    % Bode-Plot
    bode(Cpi, Cd, Cpi_ana, Cd_ana, omega_bode, opt)
    title('Cpi, Cd')
    set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

    % Legende
    if do_insert_legends
        legend('PI gemessen', 'D gemessen', 'PI analytisch', 'D analytisch')
    end
end

