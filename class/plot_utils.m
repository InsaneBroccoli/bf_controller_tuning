%==========================================================================
% PLOT_UTILS - Betaflight Controller Analysis Visualization Class
%==========================================================================
% Purpose: 
% This class provides visualization methods for analyzing quadcopter 
% flight controller performance and tuning data.
%
% Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
% Supervisor: [Michael Peter]
% Date: [25.11.2025]
%==========================================================================

classdef plot_utils

    properties
        data_flight         % Objekt of Typ flight_data
        gyro_tuning         % Objekt of Typ flight_ctrl_tuning
        analysis_flight     % Objekt of Typ flight_analyser        
        do_insert_legends    (1,1) logical = true
        linewidth            (1,1) double  = 1.2
        colorOrder = get(0, 'DefaultAxesColorOrder')
        axis_names cell = {'Roll', 'Pitch', 'Yaw'}
        pos_bode double = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
                          0.1514, 0.1100,      0.7536, 0.1917] 
    
    end
    methods
        function obj = plot_utils(data_flight, gyro_tuning, analysis_flight)
            % Konstruktor: speichert die Objekte im Plot-Helper
            obj.data_flight      = data_flight;
            obj.gyro_tuning      = gyro_tuning;
            obj.analysis_flight  = analysis_flight;
        end

  %%
        % =================================================================
        %  FIGURE 1: GYRO SIGNALS
        % =================================================================
        % Plots gyro signals for roll, pitch and yaw axes

        function obj = plot_Gyro_Signal(obj, do_insert_legends)
            df = obj.data_flight;
            figure(1);clf;
            set(gcf, 'Name', 'Flight Gyro Data');
            % --- Roll ---
            ax(1) = subplot(3,1,1); hold(ax(1),'on');
            plot(ax(1), df.time, df.data(:, df.ind.setpoint(1)));
            plot(ax(1), df.time, df.data(:, df.ind.gyroUnfilt(1)));
            plot(ax(1), df.time, df.data(:, df.ind.gyroADC(1)));
            grid(ax(1),'on');
            ylabel(ax(1),'Roll (deg/s)');
            title(ax(1),'Flight Gyro Data');
        
            if do_insert_legends
                legend(ax(1), ...
                    {'setpoint','gyro unfiltert','gyroADC'}, ...
                    'Location','best');
            end
        
            % --- Pitch ---
            ax(2) = subplot(3,1,2); hold(ax(2),'on');
            plot(ax(2), df.time, df.data(:, df.ind.setpoint(2)));
            plot(ax(2), df.time, df.data(:, df.ind.gyroUnfilt(2)));
            plot(ax(2), df.time, df.data(:, df.ind.gyroADC(2)));
            grid(ax(2),'on');
            ylabel(ax(2),'Pitch (deg/s)');
            legend(ax(2),'off');
        
            % --- Yaw ---
            ax(3) = subplot(3,1,3); hold(ax(3),'on');
            plot(ax(3), df.time, df.data(:, df.ind.setpoint(3)));
            plot(ax(3), df.time, df.data(:, df.ind.gyroUnfilt(3)));
            plot(ax(3), df.time, df.data(:, df.ind.gyroADC(3)));
            grid(ax(3),'on');
            ylabel(ax(3),'Yaw (deg/s)');
            xlabel(ax(3),'Time (s)');
            legend(ax(3),'off');
        
            % --- Axis linking & styling ---
            xmax = max(df.time);
            linkaxes(ax,'x');
            xlim(ax,[0 xmax]);
            set(findall(gcf,'Type','line'),'LineWidth',obj.linewidth);
        end

%%
        % =================================================================
        %  FIGURE 2: OVERVIEW PLOTS
        % =================================================================
        % Plots gyro signals, motor data and setpoint data

        function obj = plot_Overview(obj, do_insert_legends)
            df = obj.data_flight;
            figure(2); clf
            set(gcf, 'Name', 'Flight Overview');
                
            % ---------- Gyro Data Unfiltert ----------
            ax(1) = subplot(4,1,1);
            plot(ax(1), df.time, df.data(:, df.ind.gyroUnfilt(1:3)));
            ylabel(ax(1),'Gyro (deg/s)'); title(ax(1),'Flight Overview');
            if do_insert_legends
                legend('Roll', 'Pitch', 'Yaw', Location='best'); end
         
            % ---------- Gyro Data Axis Sum ----------
            ax(2) = subplot(4,1,2);
            plot(ax(2), df.time, df.data(:, df.ind.axisSum(1:3)));
            ylabel(ax(2),'AxisSum');
            
            % ---------- Motor Data ----------
            ax(3) = subplot(4,1,3);
            plot(ax(3), df.time, df.data(:, df.ind.motor)); 
            ylabel(ax(3),'Motor');
            if do_insert_legends
            legend('Motor 1', 'Motor 2', 'Motor 3', 'Motor 4', Location='best');
            end
        
            % ---------- Setpoint Data ----------
            ax(4) = subplot(4,1,4);
            plot(ax(4), df.time, df.data(:, df.ind.setpoint(4)));
            ylabel(ax(4),'Throttle'); xlabel(ax(4),'Time (sec)');
            grid(ax ,'on');
            
            linkaxes(ax(1:4), 'x');
             
            set(ax, 'XLim', [0, max(df.time)]);
            set(findall(gcf, 'type', 'line'), 'LineWidth', obj.linewidth);
        end


         %%
        % =================================================================
        %  FIGURE 3: Convert and evaluation time
        % =================================================================

        function obj = plot_Eval_Time(obj)
            df = obj.data_flight;
            
            figure(3); clf;
            set(gcf, 'Name', 'Evaluation Time Flight');
            delta_time_mus = diff(df.time)* 1.0e6;
            plot(df.time(1:end-1), delta_time_mus), grid on
            title(sprintf('Mean: %0.2f us, Median: %0.2f us, Std: %0.2f us\n', ...
                  mean(delta_time_mus), ...
                  median(delta_time_mus), ...
                  std(delta_time_mus)))
            xlabel('Time (sec)'), ylabel('Ts log (us)')
            xlim([0, df.time(end)])
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)            
        end

        %%
        % =================================================================
        %  FIGURE 4: GYRO SPECTRA
        % =================================================================
        % Generates spectral analysis plots for gyro data

        function obj = plot_Gyro_spectra(obj, do_insert_legends)

            fa = obj.analysis_flight;
            figure(4); clf;
            set(gcf, 'Name', 'Gyro Spectra');
            ax(1) = subplot(3, 1, 1);
            plot(ax(1), fa.freq_spectra, fa.spectra(:,1:3));
            grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)');
            set(ax(1),'YScale','log');
            title(ax(1),'Unfiltered gyro magnitude spectra');
            legend(ax(1), 'off');
            
            % ---------- Subplot 2: Filtered (ADC) gyro spectra ----------
            ax(2) = subplot(3, 1, 2);
            plot(ax(2), fa.freq_spectra, fa.spectra(:,4:6));
            grid(ax(2),'on');
            ylabel(ax(2),'Gyro (deg/s)');
            set(ax(2),'YScale','log');
            title(ax(2),'Filtered (ADC) gyro magnitude spectra');
            
            if do_insert_legends
                legend(ax(2), {'Roll','Pitch','Yaw'}, 'Location','northeast');
            end
            
            % ---------- Subplot 3: Axis sum spectra ----------
            ax(3) = subplot(3, 1, 3);
            plot(ax(3), fa.freq_spectra, fa.spectra(:,7:9));
            grid(ax(3),'on');
            ylabel(ax(3),'AxisSum');
            xlabel(ax(3),'Frequency (Hz)');
            set(ax(3),'YScale','log');
            title(ax(3),'Axis sum spectra');
            legend(ax(3), 'off');
            
            % ---------- Link axes and set Nyquist limit ----------
            linkaxes(ax, 'x');
            
            nyq = 1/(2 * fa.Ts_log);   % Nyquist frequency
            xlim([0, nyq]);
            
            % Optional uniform y-limits
            try
                ylim(ax, [1e-3 1e1]);
            catch
            end
        end

        %%
        % =================================================================
        %  FIGURE 5: GYRO SPECTRA
        % =================================================================
        % Generates spectogram plots (Thrust and Frequency)

        function obj = plot_Spectogram(obj,num_spectrograms)

            fa = obj.analysis_flight;

            figure(5);clf;
            set(gcf, 'Name', 'Gyro Spectrograms');
            sgtitle('Gyro Spectrograms')
            axes_labels = {'Roll', 'Pitch', 'Yaw'};
            c_lim = [5e-2 3e0];
            colormap('jet')

           % --- Unfiltered Spectrograms ---
            for spectrogram_nr = 1:num_spectrograms
                subplot(230 + spectrogram_nr)
                qmesh = pcolor(fa.freq_spectogram{spectrogram_nr}, ...
                               fa.throttle_all{spectrogram_nr}, ...
                               fa.spectrogram_unf{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'flat'); 
                set(qmesh, 'LineWidth', obj.linewidth);
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                xlabel('Frequency (Hz)')
                title([axes_labels{spectrogram_nr}, ' – without Filter'])   
                set(gca, 'ColorScale', 'log')                            
                clim(c_lim);                                             
                ylim([0 100])                                             
            end
        
            % --- Filtered Spectrograms ---
            for spectrogram_nr = 1:num_spectrograms
                subplot(230 + spectrogram_nr + 3)
                qmesh = pcolor(fa.freq_spectogram{spectrogram_nr}, ...
                               fa.throttle_all{spectrogram_nr}, ...
                               fa.spectrogram_fil{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'flat'); 
                set(qmesh, 'LineWidth', obj.linewidth);
                xlabel('Frequency (Hz)')
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                title([axes_labels{spectrogram_nr}, ' – with Filter'])
                colormap('jet')
                set(gca, 'ColorScale', 'log')
                clim(c_lim);
                ylim([0 100])
                set(findall(gcf, 'type', 'line'), 'LineWidth', obj.linewidth);
            end

        end
        %%
        % =================================================================
        %  FIGURE 6: BODE PLOTS
        % =================================================================
        % Bode plot (Plant and Coherence) for selected axis

        function obj = plot_Bode_Plant(obj, ind_ax)
            td = obj.gyro_tuning;
   
            figure(expand_multiple_figure_nr(6, ind_ax));clf;
            set(gcf, 'Name', ['Bode Plot - ', obj.axis_names{ind_ax}]);
            opt = bode_plot_options('dB', 'linear', 'deg', 'Hz');

            % --- Plant (Bode: magnitude + phase) ---
            ax(1) = subplot('Position', obj.pos_bode(1,:));
            bode(ax(1), td.P{ind_ax}, 'k', td.omega_bode, opt);
            title(['Plant P - ', obj.axis_names{ind_ax}]);
          
            % --- Coherence (magnitude only) ---
            ax(2) = subplot('Position', obj.pos_bode(2,:));
            opt_coh = bode_plot_options('abs', 'linear', 'deg', 'Hz');

            bodemag(ax(2), td.Coh{ind_ax}, td.omega_bode,'-k', opt_coh);
            title(''); ylabel('Coherence'); ylim([0 1]);
            linkaxes(ax, 'x'),xlim('auto'),
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

        end
%%
         % =================================================================
        %  FIGURE 7: BODE PLOTS PI & D CONTROLLER
        % =================================================================
        % Bode plot for PI and D controller (measured vs analytical)

        function obj = plot_Bode_Contr(obj, ind_ax, do_insert_legends)
            
            td = obj.gyro_tuning;
            
            opt = bode_plot_options('dB', 'linear', 'deg', 'Hz');
            
            figure(expand_multiple_figure_nr(7, ind_ax));clf;
            set(gcf, 'Name', ['Bode Controller - ', obj.axis_names{ind_ax}]);
            bode(td.Cpi{ind_ax}, td.Cpi_ana{ind_ax},td.Cd{ind_ax}, ...
                td.Cd_ana{ind_ax}, td.omega_bode, opt)              
            if do_insert_legends
                legend('PI measured','PI analytically', 'D measured', 'D analytically')
            end
            title(['Bode PI and D Controller - ', obj.axis_names{ind_ax}]);
            xlim('auto');
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

        end
        %%
        % =================================================================
        %  FIGURE 8: Gang of Four
        % =================================================================
        % Bode of Tacking, Sensitivity, Controller Effort and Compliance

        function obj = plot_Gang_of_Four(obj, do_insert_legends)
            td = obj.gyro_tuning;
            opt = bode_plot_options('dB', 'linear', 'deg', 'Hz');
            
            ind_ax = td.ind_ax;

            figure(expand_multiple_figure_nr(8, ind_ax))
            set(gcf, 'Name', ['Gang of Four - ', obj.axis_names{ind_ax}]);
            
            % --- T: Tracking ---------------
            ax(1) = subplot(2,2,1);
            bodemag(ax(1), td.CloLoAan.T, td.CloLoAanNew.T, td.T{ind_ax}, opt);
            title(ax(1), 'Tracking T');
            if do_insert_legends
                legend('actual','new','measured', 'Location','best');end

            % --- S: Sensitivity ---------
            ax(2) = subplot(2,2,2);
            bodemag(ax(2), td.CloLoAan.S, td.CloLoAanNew.S , td.omega_bode, opt); 
            title('Sensitivity S')
    
            % --- SC: Controller Effort -----------
            ax(3) = subplot(223);
            bodemag(ax(3), td.CloLoAan.SC, td.CloLoAanNew.SC, td.omega_bode, opt); 
            title('Controller Effort SC')

            % --- SP: Disturbance / Compliance -------
            ax(4) = subplot(224);
            bodemag(ax(4), td.CloLoAan.SP, td.CloLoAanNew.SP, td.omega_bode, opt); 
            title('Compliance SP')
        
            linkaxes(ax, 'x'),xlim('auto'),,
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth);
            sgtitle(['Gang of Four - ', obj.axis_names{ind_ax}]);
        end
%%
 %%
        % =================================================================
        %  FIGURE 9: Step Response
        % =================================================================
        % Step Response of choicen axis

        function obj = plot_Step_Response(obj, do_insert_legends)

            td = obj.gyro_tuning;
            
            figure(expand_multiple_figure_nr(9, td.ind_ax))
            set(gcf, 'Name', ['Step Response - ', obj.axis_names{td.ind_ax}]);
            % --- Step Response Plot ---
            ax(1) = subplot(2,1,1);
            plot(ax(1), td.step_time, td.step_resp_tra), grid on, ylabel('Gyro (deg/sec)')
            if do_insert_legends, legend('actual calculated', ...
                    'new calculated', 'measured', 'location', 'best'), end
            ylim(ax(1), 'auto');
            title(['Tracking T - ', obj.axis_names{td.ind_ax}]);
            
             % --- Compliance-Plot (SP) ---
            ax(2) = subplot(2,1,2);
            plot(ax(2), td.step_time, td.step_resp_com), grid on
            ylim(ax(2), 'auto')
            title('Compliance SP'), xlabel('Time (sec)'), ylabel('Gyro (deg/sec)')
            linkaxes(ax, 'x'), clear ax, xlim([0 0.5])
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

        end

    end
end