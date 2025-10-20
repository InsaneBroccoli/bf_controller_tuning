function fig_nr = expand_multiple_figure_nr(ini_fig_nr, multp_fig_nr)
%EXPAND_MULTIPLE_FIGURE_NR  Expand a figure number into repeated digits
%   fig_nr = expand_multiple_figure_nr(ini_fig_nr, multp_fig_nr)
%
% PURPOSE
%   - Create a new figure number by repeating the digits of the initial number
%   - Used to generate distinct but related figure identifiers
%
% INPUTS
%   - ini_fig_nr   scalar base figure number
%   - multp_fig_nr scalar multiplier, number of repetitions
%
% OUTPUTS
%   - fig_nr scalar expanded figure number
%
% METHOD
%   1) Convert ini_fig_nr to string
%   2) Concatenate the string multp_fig_nr times
%   3) Convert back to number
%
% NOTES
%   - If multp_fig_nr == 1, the function simply returns ini_fig_nr
%   - Example: ini_fig_nr = 12, multp_fig_nr = 3 → fig_nr = 121212

    if multp_fig_nr == 1
        fig_nr = ini_fig_nr;
        return
    else
        fig_nr_str = num2str(ini_fig_nr);
        for i = 1:multp_fig_nr-1
            fig_nr_str = [fig_nr_str, num2str(ini_fig_nr)]; %#ok
        end
    end
    fig_nr = str2double(fig_nr_str);

end
