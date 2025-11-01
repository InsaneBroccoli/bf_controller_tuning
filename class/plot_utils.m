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
        opt
        colorOrder = get(0, 'DefaultAxesColorOrder')
        second_flight (1,1) logical = false
        ind_ax
        pos_bode double = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
        0.1514, 0.1100,      0.7536, 0.1917] 

    end


    methods
  %%
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
                        'Location','northeastoutside');
         
                else
                    legend(ax(1), ...
                        {'setpoint F1','gyro F1','gyroADC F1'}, ...
                        'Location','northeastoutside','Orientation','horizontal');
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

%%
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
                        'Location','northeastoutside');
                else
                    legend(ax(1), {'F1 Roll','F1 Pitch','F1 Yaw'}, 'Location','northeastoutside');
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
%%
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
            legend('Roll', 'Pitch', 'Yaw', Location='northwest');
        
            ax(2) = subplot(4, ncol, pos(2,1));
            plot(ax(2), flight1.time, flight1.axisSumData); grid(ax(2),'on');
            ylabel(ax(2),'AxisSum');
        
            ax(3) = subplot(4, ncol, pos(3,1));
            plot(ax(3), flight1.time, flight1.motorData); grid(ax(3),'on');
            ylabel(ax(3),'Motor');
            legend('Motor 1', 'Motor 2', 'Motor 3', 'Motor 4', Location='north');
        
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

%%
        % ===============================================================
        %  FIGURE 4: B PLOTS
        % ===============================================================

        function plotBode(obj, flight1, varargin)
            % plotBode - Bode plot (Plant and Coherence) for selected axis
            % Minimal edits: keep all global bode options in MAIN.
            % Use local options only for the coherence subplot (auto Y-limits).
            % Always use BODEMAG for coherence on both flights.
        
            if obj.second_flight
                flight2 = varargin{1};
            end
        
            switch obj.ind_ax
                case 1
                    figure(4)
                    ax(1) = subplot('Position', obj.pos_bode(1,:));
        
                    % --- Plant (Bode: magnitude + phase) ---
                    bode(ax(1), flight1.transfData, 'k', flight1.transfOmega, obj.opt), title('Plant P Roll')
                    if obj.second_flight
                        hold on
                        bode(ax(1), flight2.transfData, 'r', flight2.transfOmega, obj.opt)
                    end
                    hold off, grid on
                    legend('Flight 1', 'Flight 2', 'Location','southwest')
        
                    % --- Coherence (magnitude only) ---
                    ax(2) = subplot('Position', obj.pos_bode(2,:));
                    
                    bodemag(ax(2), flight1.transfCoher, 'k', flight1.transfOmega, obj.opt),
                    title(''), ylabel('Coherence')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.transfCoher, 'r', flight2.transfOmega, obj.opt)
                    end
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                case 2
                    figure(44)
                    ax(1) = subplot('Position', obj.pos_bode(1,:));
        
                    % --- Plant (Bode: magnitude + phase) ---
                    bode(ax(1), flight1.transfData, 'k', flight1.transfOmega, obj.opt), title('Plant P Pitch')
                    if obj.second_flight
                        hold on
                        bode(ax(1), flight2.transfData, 'r', flight2.transfOmega, obj.opt)
                    end
                    hold off, grid on
                    legend('Flight 1', 'Flight 2', 'Location','southwest')
        
                    % --- Coherence (magnitude only) ---
                    ax(2) = subplot('Position', obj.pos_bode(2,:));
                    bodemag(ax(2), flight1.transfCoher, 'k', flight1.transfOmega, obj.opt), title(''), ylabel('Coherence')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.transfCoher, 'r', flight2.transfOmega, obj.opt)
                    end
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                case 3
                    figure(4)
                    ax(1) = subplot('Position', obj.pos_bode(1,:));
        
                    % --- Plant (Bode: magnitude + phase) ---
                    bode(ax(1), flight1.transfData, 'k', flight1.transfOmega, obj.opt), title('Plant P Yaw')
                    if obj.second_flight
                        hold on
                        bode(ax(1), flight2.transfData, 'r', flight2.transfOmega, obj.opt)
                    end
                    hold off, grid on
                    legend('Flight 1', 'Flight 2', 'Location','southwest')
        
                    % --- Coherence (magnitude only) ---
                    ax(2) = subplot('Position', obj.pos_bode(2,:));
                    bodemag(ax(2), flight1.transfCoher, 'k', flight1.transfOmega, obj.opt), title(''), ylabel('Coherence')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.transfCoher, 'r', flight2.transfOmega, obj.opt)
                    end
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
            end
        end

    %%
        % ===============================================================
        %  FIGURE 5: CONTROLLER BODE PLOTS
        % ===============================================================

        function plotCPIDBode(obj, flight1, varargin)
            % plotCPIDBode - Bode plot for PI and D controller (measured vs analytical)
            % Minimal edits: keep global bode options from MAIN, do not override here.
        
            if obj.second_flight
                flight2 = varargin{1};
            end
        
            switch obj.ind_ax
                case 1
                    figure(5)
        
                    % ---------- PI ----------
                    subplot(1,2,1)
                    bode(flight1.transfCpi, flight1.transfCpiAna, flight1.transfOmega, obj.opt), title('Controller PI Roll')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCpi, flight2.transfCpiAna, flight2.transfOmega, obj.opt)
                        
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1', ...
                                   'PI gemessen F2','PI analytisch F2')
                        end
                    else
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                    % ---------- D ----------
                    subplot(1,2,2)
                    bode(flight1.transfCD, flight1.transfCDAna, flight1.transfOmega, obj.opt), title('Controller D Roll')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCD, flight2.transfCDAna, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1', ...
                                   'D gemessen F2','D analytisch F2', ...
                                   'Location','northwest')
                        end
                    else
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                case 2
                    figure(55)
        
                    % ---------- PI ----------
                    subplot(1,2,1)
                    bode(flight1.transfCpi, flight1.transfCpiAna, flight1.transfOmega, obj.opt), title('Controller PI Pitch')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCpi, flight2.transfCpiAna, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1', ...
                                   'PI gemessen F2','PI analytisch F2')
                        end
                    else
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                    % ---------- D ----------
                    subplot(1,2,2)
                    bode(flight1.transfCD, flight1.transfCDAna, flight1.transfOmega, obj.opt), title('Controller D Pitch')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCD, flight2.transfCDAna, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1', ...
                                   'D gemessen F2','D analytisch F2', ...
                                   'Location','northwest')
                        end
                    else
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                case 3
                    figure(5)
        
                    % ---------- PI ----------
                    subplot(1,2,1)
                    bode(flight1.transfCpi, flight1.transfCpiAna, flight1.transfOmega, obj.opt), title('Controller PI Yaw')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCpi, flight2.transfCpiAna, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1', ...
                                   'PI gemessen F2','PI analytisch F2')
                        end
                    else
                        if obj.do_insert_legends
                            legend('PI gemessen F1','PI analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
        
                    % ---------- D ----------
                    subplot(1,2,2)
                    bode(flight1.transfCD, flight1.transfCDAna, flight1.transfOmega, obj.opt), title('Controller D Yaw')
                    hold on
                    if obj.second_flight
                        bode(flight2.transfCD, flight2.transfCDAna, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1', ...
                                   'D gemessen F2','D analytisch F2', ...
                                   'Location','northwest')
                        end
                    else
                        if obj.do_insert_legends
                            legend('D gemessen F1','D analytisch F1')
                        end
                    end
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
            end
        end


    %%
        % ===============================================================
        %  FIGURE 6: Gang of Four
        % ===============================================================

        function plotGangofFour(obj, flight1, varargin)
          
        if obj.second_flight
            flight2 = varargin{1};
        end
            switch obj.ind_ax
                %============ ROLL ==========
                case 1
                    figure(6)
                    ax(1) = subplot(2,2,1);
                    bodemag(ax(1), flight1.CloLoAan.T , flight1.CloLoAanNew.T , flight1.transfT, flight1.transfOmega, obj.opt), title('Tracking T ')
                    if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bodemag(ax(1), flight2.CloLoAan.T , flight2.CloLoAanNew.T , flight2.transfT, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', 'actual F2', 'new F2','location', 'best'), end
                    end
                    ax(2) = subplot(2,2,2);
                    bodemag(ax(2), flight1.CloLoAan.S, flight1.CloLoAanNew.S , flight1.transfOmega, obj.opt), title('Sensitivity S')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.CloLoAan.S, flight2.CloLoAanNew.S , flight2.transfOmega, obj.opt)
                    end
                    ax(3) = subplot(223);
                    bodemag(ax(3), flight1.CloLoAan.SC, flight1.CloLoAanNew.SC, flight1.transfOmega, obj.opt), title('Controller Effort SC')
                    if obj.second_flight
                        hold on
                        bodemag(ax(3), flight2.CloLoAan.SC, flight2.CloLoAanNew.SC, flight2.transfOmega, obj.opt)
                    end
                    ax(4) = subplot(224);
                    bodemag(ax(4), flight1.CloLoAan.SP, flight1.CloLoAanNew.SP, flight1.transfOmega, obj.opt), title('Compliance SP')
                    if obj.second_flight
                        hold on
                        bodemag(ax(4), flight2.CloLoAan.SP, flight2.CloLoAanNew.SP, flight2.transfOmega, obj.opt)
                    end
        
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
                    sgtitle('Roll')

                    %====== PITCH =======
                case 2
                    figure(66)
                    ax(1) = subplot(2,2,1);
                    bodemag(ax(1), flight1.CloLoAan.T , flight1.CloLoAanNew.T , flight1.transfT, flight1.transfOmega, obj.opt), title('Tracking T')
                    if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bodemag(ax(1), flight2.CloLoAan.T , flight2.CloLoAanNew.T , flight2.transfT, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', 'actual F2', 'new F2','location', 'best'), end
                    end
                    ax(2) = subplot(2,2,2);
                    bodemag(ax(2), flight1.CloLoAan.S, flight1.CloLoAanNew.S , flight1.transfOmega, obj.opt), title('Sensitivity S')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.CloLoAan.S, flight2.CloLoAanNew.S , flight2.transfOmega, obj.opt)
                    end
                    ax(3) = subplot(223);
                    bodemag(ax(3), flight1.CloLoAan.SC, flight1.CloLoAanNew.SC, flight1.transfOmega, obj.opt), title('Controller Effort SC')
                    if obj.second_flight
                        hold on
                        bodemag(ax(3), flight2.CloLoAan.SC, flight2.CloLoAanNew.SC, flight2.transfOmega, obj.opt)
                    end
                    ax(4) = subplot(224);
                    bodemag(ax(4), flight1.CloLoAan.SP, flight1.CloLoAanNew.SP, flight1.transfOmega, obj.opt), title('Compliance SP')
                    if obj.second_flight
                        hold on
                        bodemag(ax(4), flight2.CloLoAan.SP, flight2.CloLoAanNew.SP, flight2.transfOmega, obj.opt)
                    end
        
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
                    sgtitle('Pitch')

                    %========== YAW =========
                case 3
                    figure(666)
                    ax(1) = subplot(2,2,1);
                    bodemag(ax(1), flight1.CloLoAan.T , flight1.CloLoAanNew.T , flight1.transfT, flight1.transfOmega, obj.opt), title('Tracking T')
                    if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bodemag(ax(1), flight2.CloLoAan.T , flight2.CloLoAanNew.T , flight2.transfT, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', 'actual F2', 'new F2','location', 'best'), end
                    end
                    ax(2) = subplot(2,2,2);
                    bodemag(ax(2), flight1.CloLoAan.S, flight1.CloLoAanNew.S , flight1.transfOmega, obj.opt), title('Sensitivity S')
                    if obj.second_flight
                        hold on
                        bodemag(ax(2), flight2.CloLoAan.S, flight2.CloLoAanNew.S , flight2.transfOmega, obj.opt)
                    end
                    ax(3) = subplot(223);
                    bodemag(ax(3), flight1.CloLoAan.SC, flight1.CloLoAanNew.SC, flight1.transfOmega, obj.opt), title('Controller Effort SC')
                    if obj.second_flight
                        hold on
                        bodemag(ax(3), flight2.CloLoAan.SC, flight2.CloLoAanNew.SC, flight2.transfOmega, obj.opt)
                    end
                    ax(4) = subplot(224);
                    bodemag(ax(4), flight1.CloLoAan.SP, flight1.CloLoAanNew.SP, flight1.transfOmega, obj.opt), title('Compliance SP')
                    if obj.second_flight
                        hold on
                        bodemag(ax(4), flight2.CloLoAan.SP, flight2.CloLoAanNew.SP, flight2.transfOmega, obj.opt)
                    end
        
                    linkaxes(ax, 'x'), clear ax
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
                    sgtitle('Yaw')
            end
        end


            %%
        % ===============================================================
        %  FIGURE 7: Step Response
        % ===============================================================

        function plotStepResp(obj, flight1, varargin)
            if obj.second_flight
                flight2 = varargin{1};
            end
            switch obj.ind_ax
                case 1
            % ====== ROLL ======
                figure(7)
                ax(1) = subplot(2,1,1);
                plot(ax(1), flight1.step_resp_tim, flight1.step_resp_tra), grid on, ylabel('Gyro (deg/sec)')
                if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                if obj.second_flight
                    hold on
                    plot(ax(1), flight2.step_resp_tim, flight2.step_resp_tra)
                    if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                            'actual F2', 'new F2','location', 'best'), end
                end
                title('Tracking T Roll')
                ylim([0 1.3])
    
                ax(2) = subplot(2,1,2);
                plot(ax(2), flight1.step_resp_tim, flight1.step_resp_com), grid on
                if obj.second_flight
                    hold on
                    plot(ax(2), flight2.step_resp_tim, flight2.step_resp_tra)
                end
                title('Compliance SP Roll'), xlabel('Time (sec)'), ylabel('Gyro (deg/sec)')
                ylim([-0.2 1.1])
                linkaxes(ax, 'x'), clear ax, xlim([0 0.5])
                set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

                % ====== PITCH ======
                case 2
                    figure(77)
                    ax(1) = subplot(2,1,1);
                    plot(ax(1), flight1.step_resp_tim, flight1.step_resp_tra), grid on, ylabel('Gyro (deg/sec)')
                    if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        plot(ax(1), flight2.step_resp_tim, flight2.step_resp_tra)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                                'actual F2', 'new F2','location', 'best'), end
                    end
                    title('Tracking T Pitch')
                    ylim([0 1.3])
        
                    ax(2) = subplot(2,1,2);
                    plot(ax(2), flight1.step_resp_tim, flight1.step_resp_com), grid on
                    if obj.second_flight
                        hold on
                        plot(ax(2), flight2.step_resp_tim, flight2.step_resp_tra)
                    end
                    title('Compliance SP Pitch'), xlabel('Time (sec)'), ylabel('Gyro (deg/sec)')
                    ylim([-0.2 1.1])
                    linkaxes(ax, 'x'), clear ax, xlim([0 0.5])
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

                % ====== YAW ======
                case 3

                    figure(777)
                    ax(1) = subplot(2,1,1);
                    plot(ax(1), flight1.step_resp_tim, flight1.step_resp_tra), grid on, ylabel('Gyro (deg/sec)')
                    if obj.do_insert_legends, legend('actual', 'new', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        plot(ax(1), flight2.step_resp_tim, flight2.step_resp_tra)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                                'actual F2', 'new F2','location', 'best'), end
                    end
                    title('Tracking T Pitch')
                    ylim([0 1.3])
        
                    ax(2) = subplot(2,1,2);
                    plot(ax(2), flight1.step_resp_tim, flight1.step_resp_com), grid on
                    if obj.second_flight
                        hold on
                        plot(ax(2), flight2.step_resp_tim, flight2.step_resp_tra)
                    end
                    title('Compliance SP Pitch'), xlabel('Time (sec)'), ylabel('Gyro (deg/sec)')
                    ylim([-0.2 1.1])
                    linkaxes(ax, 'x'), clear ax, xlim([0 0.5])
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

            end
        end
    %%
        % ===============================================================
        %  FIGURE 8: Bode Controller
        % ===============================================================

        function plotController(obj, flight1, varargin)
            if obj.second_flight
                flight2 = varargin{1};
            end
            switch obj.ind_ax
                %======== ROLL ========
                case 1
                    figure(8)
                    bode(flight1.CloLoAan.C, flight1.CloLoAanNew.C, flight1.transfOmega, obj.opt)
                    if obj.do_insert_legends, legend('actual F1', 'new F1', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bode(flight2.CloLoAan.C, flight2.CloLoAanNew.C, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                                'actual F2', 'new F2', 'location', 'best'), end
                    end
                   
                    title('Controller C Roll')
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
                
                %======== PITCH ========
                case 2
                    figure(88)
                    bode(flight1.CloLoAan.C, flight1.CloLoAanNew.C, flight1.transfOmega, obj.opt)
                    if obj.do_insert_legends, legend('actual F1', 'new F1', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bode(flight2.CloLoAan.C, flight2.CloLoAanNew.C, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                                'actual F2', 'new F2', 'location', 'best'), end
                    end

                    title('Controller C Pitch')
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
                %======== YAW ========
                case 3
                    figure(88)
                    bode(flight1.CloLoAan.C, flight1.CloLoAanNew.C, flight1.transfOmega, obj.opt)
                    if obj.do_insert_legends, legend('actual F1', 'new F1', 'location', 'best'), end
                    if obj.second_flight
                        hold on
                        bode(flight2.CloLoAan.C, flight2.CloLoAanNew.C, flight2.transfOmega, obj.opt)
                        if obj.do_insert_legends, legend('actual F1', 'new F1', ...
                                'actual F2', 'new F2', 'location', 'best'), end
                    end
                    
                    title('Controller C Yaw')
                    set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
            end

        end


         %%
        % ===============================================================
        %  FIGURE 98/99: Convert and evaluate time
        % ===============================================================

        function plotevaltime(obj, flight1, varargin)
            if obj.second_flight
                flight2 = varargin{1};
            end

            delta_time_mus = diff(flight1.time)* 1.0e6;
            figure(98); clf
            plot(flight1.time(1:end-1), delta_time_mus), grid on
            title(sprintf('Flight 1: Mean: %0.2f mus, Median: %0.2f mus, Std: %0.2f mus\n', ...
                  mean(delta_time_mus), ...
                  median(delta_time_mus), ...
                  std(delta_time_mus)))
            xlabel('Time (sec)'), ylabel('Ts log (mus)')
            xlim([0, flight1.time(end)])
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)

            if obj.second_flight
                delta_time_mus = diff(flight2.time)* 1.0e6;
            figure(99); clf
            plot(flight2.time(1:end-1), delta_time_mus), grid on
            title(sprintf('Flight 2: Mean: %0.2f mus, Median: %0.2f mus, Std: %0.2f mus\n', ...
                  mean(delta_time_mus), ...
                  median(delta_time_mus), ...
                  std(delta_time_mus)))
            xlabel('Time (sec)'), ylabel('Ts log (mus)')
            xlim([0, flight2.time(end)])
            set(findall(gcf, 'type', 'line'), 'linewidth', obj.linewidth)
            end
            
        end

         %%
        % ===============================================================
        %  FIGURE 22/23: Spectogram
        % ===============================================================

        function plotspectogram(obj, flight1, varargin)

            figure(22); clf
            sgtitle('Gyro Spectrograms Flight 1')
            axes_labels = {'Roll', 'Pitch', 'Yaw'};
        
            c_lim = [5e-2 3e0];
        
            colormap('jet')
            
           for spectrogram_nr = 1:flight1.num_spectrograms
                subplot(230 + spectrogram_nr)
                qmesh = pcolor(flight1.SpecfreqUnf{spectrogram_nr}, ...
                               flight1.SpecthroUnf{spectrogram_nr}, ...
                               flight1.SpecAmplUnf{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                xlabel('Frequency (Hz)')
                title([axes_labels{spectrogram_nr}, ' – ohne Filter'])   
                set(gca, 'ColorScale', 'log')                            
                clim(c_lim);                                             
                ylim([0 100])                                             
            end
        
            % --- unten: mit Filter (unverändert) ---
            for spectrogram_nr = 1:flight1.num_spectrograms
                subplot(230 + spectrogram_nr + 3)
                qmesh = pcolor(flight1.SpecfreqFil{spectrogram_nr}, ...
                               flight1.SpecthroFil{spectrogram_nr}, ...
                               flight1.SpecAmplFil{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
                xlabel('Frequency (Hz)')
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                title([axes_labels{spectrogram_nr}, ' – mit Filter'])
                colormap('jet')
                set(gca, 'ColorScale', 'log')
                clim(c_lim);
                ylim([0 100])
            end

            % ======== Spectogram Flight 2 ========

            if obj.second_flight
                flight2 = varargin{1};
            end
            figure(23); clf
            sgtitle('Gyro Spectrograms Flight 2')
            axes_labels = {'Roll', 'Pitch', 'Yaw'};
        
            c_lim = [5e-2 3e0];
        
            colormap('jet')
            
           for spectrogram_nr = 1:flight2.num_spectrograms
                subplot(230 + spectrogram_nr)
                qmesh = pcolor(flight2.SpecfreqUnf{spectrogram_nr}, ...
                               flight2.SpecthroUnf{spectrogram_nr}, ...
                               flight2.SpecAmplUnf{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                xlabel('Frequency (Hz)')
                title([axes_labels{spectrogram_nr}, ' – ohne Filter'])   
                set(gca, 'ColorScale', 'log')                            
                clim(c_lim);                                             
                ylim([0 100])                                             
            end
        
            % --- unten: mit Filter (unverändert) ---
            for spectrogram_nr = 1:flight2.num_spectrograms
                subplot(230 + spectrogram_nr + 3)
                qmesh = pcolor(flight2.SpecfreqFil{spectrogram_nr}, ...
                               flight2.SpecthroFil{spectrogram_nr}, ...
                               flight2.SpecAmplFil{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
                xlabel('Frequency (Hz)')
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                title([axes_labels{spectrogram_nr}, ' – mit Filter'])
                colormap('jet')
                set(gca, 'ColorScale', 'log')
                clim(c_lim);
                ylim([0 100])
            end

        end

    end
end