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
        function plotGyroSpectra(obj, flight1)
            figure(2); clf

            ax(1) = subplot(2,1,1);
            plot(ax(1), flight1.Specfreq, flight1.unfgyroSpec);
            grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)');
            set(ax(1),'YScale','log');
            title(ax(1),'Magnitude Spectra');

            if obj.do_insert_legends
                legend(ax(1), ...
                    {'gyro Roll','gyro Pitch','gyro Yaw','gyroADC Roll','gyroADC Pitch','gyroADC Yaw'}, ...
                    'Location','best');
            end

            ax(2) = subplot(2,1,2);
            plot(ax(2), flight1.Specfreq, flight1.axisSumSpec);
            grid(ax(2),'on');
            ylabel(ax(2),'AxisSum');
            xlabel(ax(2),'Frequency (Hz)');
            set(ax(2),'YScale','log');

            if obj.do_insert_legends
                legend(ax(2), {'axisSum Roll','axisSum Pitch','axisSum Yaw'}, 'Location','best');
            end

            linkaxes(ax, 'x');
            xlim(ax, [0, 1/(2*flight1.Ts_log)]);
            ylim(ax(1), [1e-3 1e1]);
            ylim(ax(2), [1e-3 1e1]);
            set(findall(gcf, 'type', 'line'), 'LineWidth', obj.linewidth);
        end


        % ===============================================================
        %  FIGURE 3: OVERVIEW PLOTS
        % ===============================================================
        function plotOverview(obj, flight1)
            figure(3); clf

            ax(1) = subplot(4,1,1);
            plot(ax(1), flight1.time, flight1.unfgyroData);
            grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)');

            ax(2) = subplot(4,1,2);
            plot(ax(2), flight1.time, flight1.axisSumData);
            grid(ax(2),'on');
            ylabel(ax(2),'AxisSum');

            ax(3) = subplot(4,1,3);
            plot(ax(3), flight1.time, flight1.motorData);
            grid(ax(3),'on');
            ylabel(ax(3),'Motor');

            ax(4) = subplot(4,1,4);
            plot(ax(4), flight1.time, flight1.setpoint(:,4));
            grid(ax(4),'on');
            ylabel(ax(4),'Throttle');
            xlabel(ax(4),'Time (sec)');

            linkaxes(ax,'x');
            xlim(ax,[0, flight1.time(end)]);
            set(findall(gcf,'type','line'),'linewidth',obj.linewidth);
        end
    end
end
