
function show_controller_analysis(CL_ana, CL_ana_new, T, omega_bode, do_insert_legends)

% Closed-loop bode plots (gang of four)
figure(expand_multiple_figure_nr(6, multp_fig_nr))
ax(1) = subplot(221);
opt.YLim = {[1e-3 1e1], [-180 180]}; opt.MagScale = 'log';
bodemag(ax(1), CL_ana.T , CL_ana_new.T , T, omega_bode, opt), title('Tracking T')
if do_insert_legends, legend('actual', 'new', 'location', 'best'), end
ax(2) = subplot(222);
bodemag(ax(2), CL_ana.S , CL_ana_new.S , omega_bode, opt), title('Sensitivity S')
ax(3) = subplot(223);
opt.YLim = {[1e-2 1e2], [-180 180]};
bodemag(ax(3), CL_ana.SC, CL_ana_new.SC, omega_bode, opt), title('Controller Effort SC')
ax(4) = subplot(224);
opt.YLim = {[1e-3 1e1], [-180 180]};
bodemag(ax(4), CL_ana.SP, CL_ana_new.SP, omega_bode, opt), title('Compliance SP')
linkaxes(ax, 'x'), clear ax
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

end
