classdef plot_utils < handle
    % PLOT_UTILS - simple helper class for all plotting functions
    % -----------------------------------------------------------------
    % Keeps global plot settings consistent across all figures
    % and provides reusable plotting methods.
    % -----------------------------------------------------------------

    properties
        do_compensate_iterm  (1,1) logical = true
        do_show_spec_figures (1,1) logical = true
        do_insert_legends    (1,1) logical = true
        linewidth            (1,1) double  = 1.2
        opt = bodeoptions('cstprefs')
        colorOrder = get(0, 'DefaultAxesColorOrder')
        second_flight (1,1) logical = false
        ind_ax
        pos_bode double = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
        0.1514, 0.1100,      0.7536, 0.1917] 

    end


    methods
        % ===============================================================
        %  FIGURE 1: GYRO SIGNALS
        % ===============================================================
        function plotGyroSignals(obj, flight1, varargin)
            % plotGyroSignals(flight1)
            % plotGyroSignals(flight1, flight2)
            % Creates Roll, Pitch, Yaw subplots with optional second flight
            if obj.second_flight
                flight2 = varargin{1};
            end

            figure(1); clf

            % --- Roll ---
            ax(1) = subplot(3,1,1); hold(ax(1),'on');
            plot(ax(1), flight1.time, flight1.setpoint(:,1),   'DisplayName','setpoint F1');
            plot(ax(1), flight1.time, flight1.unfgyroData(:,1),'DisplayName','gyro F1');
            plot(ax(1), flight1.time, flight1.gyroData(:,1),   'DisplayName','gyroADC F1');

            if obj.second_flight
                plot(ax(1), flight2.time, flight2.setpoint(:,1),    'DisplayName','setpoint F2');
                plot(ax(1), flight2.time, flight2.unfgyroData(:,1), 'DisplayName','gyro F2');
                plot(ax(1), flight2.time, flight2.gyroData(:,1),    'DisplayName','gyroADC F2');
            end

            grid(ax(1),'on');
            ylabel(ax(1),'Roll (deg/s)');
            title(ax(1),'Gyro Signals');

            if obj.do_insert_legends
                if obj.second_flight
                    legend(ax(1), ...
                        {'setpoint F1','gyro F1','gyroADC F1', ...
                         'setpoint F2','gyro F2','gyroADC F2'}, ...
                         'Location','northoutside','Orientation','horizontal');
                else
                    legend(ax(1), ...
                        {'setpoint F1','gyro F1','gyroADC F1'}, ...
                        'Location','northoutside','Orientation','horizontal');
                end
            end

            % --- Pitch ---
            ax(2) = subplot(3,1,2); hold(ax(2),'on');
            plot(ax(2), flight1.time, flight1.setpoint(:,2));
            plot(ax(2), flight1.time, flight1.unfgyroData(:,2));
            plot(ax(2), flight1.time, flight1.gyroData(:,2));
            if obj.second_flight
                plot(ax(2), flight2.time, flight2.setpoint(:,2));
                plot(ax(2), flight2.time, flight2.unfgyroData(:,2));
                plot(ax(2), flight2.time, flight2.gyroData(:,2));
            end
            grid(ax(2),'on');
            ylabel(ax(2),'Pitch (deg/s)');
            legend(ax(2),'off');

            % --- Yaw ---
            ax(3) = subplot(3,1,3); hold(ax(3),'on');
            plot(ax(3), flight1.time, flight1.setpoint(:,3));
            plot(ax(3), flight1.time, flight1.unfgyroData(:,3));
            plot(ax(3), flight1.time, flight1.gyroData(:,3));
            if obj.second_flight
                plot(ax(3), flight2.time, flight2.setpoint(:,3));
                plot(ax(3), flight2.time, flight2.unfgyroData(:,3));
                plot(ax(3), flight2.time, flight2.gyroData(:,3));
            end
            grid(ax(3),'on');
            ylabel(ax(3),'Yaw (deg/s)');
            xlabel(ax(3),'Time (s)');
            legend(ax(3),'off');

            % --- Axis linking & styling ---
            if obj.second_flight
                xmax = max([max(flight1.time), max(flight2.time)]);
            else
                xmax = max(flight1.time);
            end
            linkaxes(ax,'x');
            xlim(ax,[0 xmax]);
            set(findall(gcf,'Type','line'),'LineWidth',obj.linewidth);
        end


        % ===============================================================
        %  FIGURE 2: GYRO SPECTRA
        % ===============================================================
        function plotGyroSpectra(obj, flight1, varargin)
            figure(2); clf
        
            % Detect optional second flight
            if obj.second_flight
                flight2 = varargin{1};
            end
        
            % Preallocate 3 axes (3 rows × 1 col)
            ax = gobjects(1,3);
        
            % ---------- Subplot 1: Unfiltered spectra ----------
            ax(1) = subplot(3,1,1);
            plot(ax(1), flight1.Specfreq, flight1.unfgyroSpec); % F1
            hold(ax(1),'on');
            if obj.second_flight
                plot(ax(1), flight2.Specfreq, flight2.unfgyroSpec); % F2
            end
            grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)');
            set(ax(1),'YScale','log');
            title(ax(1),'Unfiltered gyro magnitude spectra');
        
            if obj.do_insert_legends
                if obj.second_flight
                    legend(ax(1), ...
                        {'F1 Roll','F1 Pitch','F1 Yaw', ...
                         'F2 Roll','F2 Pitch','F2 Yaw'}, ...
                        'Location','northeast');
                else
                    legend(ax(1), {'F1 Roll','F1 Pitch','F1 Yaw'}, 'Location','northeast');
                end
            end
        
            % ---------- Subplot 2: Filtered (ADC) spectra ----------
            ax(2) = subplot(3,1,2);
            plot(ax(2), flight1.Specfreq, flight1.adcgyroSpec); % F1
            hold(ax(2),'on');
            if obj.second_flight
                plot(ax(2), flight2.Specfreq, flight2.adcgyroSpec); % F2
            end
            grid(ax(2),'on');
            ylabel(ax(2),'Gyro (deg/s)');
            set(ax(2),'YScale','log');
            title(ax(2),'Filtered (ADC) gyro magnitude spectra');
            legend off;
          
            % ---------- Subplot 3: AxisSum spectra ----------
            ax(3) = subplot(3,1,3);
            plot(ax(3), flight1.Specfreq, flight1.axisSumSpec); % F1
            hold(ax(3),'on');
            if obj.second_flight
                plot(ax(3), flight2.Specfreq, flight2.axisSumSpec); % F2
            end
            grid(ax(3),'on');
            ylabel(ax(3),'AxisSum');
            xlabel(ax(3),'Frequency (Hz)');
            set(ax(3),'YScale','log');
            title(ax(3),'Axis sum spectra');
            legend off;

            % Link x-axes and set limits based on (smallest) Nyquist
            linkaxes(ax,'x');
            nyq1 = 1/(2*flight1.Ts_log);
            if obj.second_flight && isfield(flight2,'Ts_log')
                nyq2 = 1/(2*flight2.Ts_log);
            else
                nyq2 = nyq1;
            end
            xlim(ax,[0, min(nyq1, nyq2)]);
        
            % Optional y-limits (keep your previous styling)
            try
                ylim(ax(1), [1e-3 1e1]);
                ylim(ax(2), [1e-3 1e1]);
                ylim(ax(3), [1e-3 1e1]);
            catch
                % ignore if ranges are incompatible
            end
        
            % Apply global line width
            set(findall(gcf,'type','line'), 'LineWidth', obj.linewidth);
        end

        % ===============================================================
        %  FIGURE 3: OVERVIEW PLOTS
        % ===============================================================
        function plotOverview(obj, flight1, varargin)
            figure(3); clf
        
            ncol = 1;
            if obj.second_flight
                ncol = 2;
                flight2 = varargin{1};
            end
        
            % helper function for subplot position
            pos = @(row, colIdx) colIdx + (row-1)*ncol;   % row = 1..4, colIdx = 1 or 2
        
            % ---------- Flight 1 (left column) ----------
            ax(1) = subplot(4, ncol, pos(1,1));
            plot(ax(1), flight1.time, flight1.unfgyroData); grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)'); title(ax(1),'Flight 1');
        
            ax(2) = subplot(4, ncol, pos(2,1));
            plot(ax(2), flight1.time, flight1.axisSumData); grid(ax(2),'on');
            ylabel(ax(2),'AxisSum');
        
            ax(3) = subplot(4, ncol, pos(3,1));
            plot(ax(3), flight1.time, flight1.motorData); grid(ax(3),'on');
            ylabel(ax(3),'Motor');
        
            ax(4) = subplot(4, ncol, pos(4,1));
            plot(ax(4), flight1.time, flight1.setpoint(:,4)); grid(ax(4),'on');
            ylabel(ax(4),'Throttle'); xlabel(ax(4),'Time (sec)');
        
            % ---------- Flight 2 (right column, optional) ----------
            if obj.second_flight
                ax(5) = subplot(4, ncol, pos(1,2));
                plot(ax(5), flight2.time, flight2.unfgyroData); grid(ax(5),'on');
                ylabel(ax(5),'Gyro (deg/s)'); title(ax(5),'Flight 2');
        
                ax(6) = subplot(4, ncol, pos(2,2));
                plot(ax(6), flight2.time, flight2.axisSumData); grid(ax(6),'on');
                ylabel(ax(6),'AxisSum');
        
                ax(7) = subplot(4, ncol, pos(3,2));
                plot(ax(7), flight2.time, flight2.motorData); grid(ax(7),'on');
                ylabel(ax(7),'Motor');
        
                ax(8) = subplot(4, ncol, pos(4,2));
                plot(ax(8), flight2.time, flight2.setpoint(:,4)); grid(ax(8),'on');
                ylabel(ax(8),'Throttle'); xlabel(ax(8),'Time (sec)');
            end
        
            % ---------- Axis linking and styling ----------
            if obj.second_flight
                linkaxes(ax(1:8), 'x');
            else
                linkaxes(ax(1:4), 'x');
            end
            set(ax, 'XLim', [0, max(flight1.time)]);
            set(findall(gcf, 'type', 'line'), 'LineWidth', obj.linewidth);
        end

    end
end