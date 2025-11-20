classdef flight_data
    properties
        % File handling
        file_path        % Path to the flight log file (.bbl.csv)
        
        % Raw flight data
        time             % Time vector [s]
        data             % Flgiht data
        ind              % Columns
        para             % Parameters old
        
        % Calculated data
        Ts_log           % Logging Time [s]
        Ts_cntr          

        % Default values
        linewidth            (1,1) double  = 1.2
        
    end

    methods

        function obj = flight_data(file_path)
            obj.file_path = file_path;
        end

         function obj = get_data(obj)
            % --- Load and Process Flight Log Data ---
            [obj.para, Nheader, obj.ind, ind_cntr] = extract_header_information(obj.file_path);
            
            % Load data from CSV or cached MAT file for faster processing
            tic
            try
               load_data = load([obj.file_path(1:end-8), '.mat']);
               obj.data = load_data.data;
            catch exception
               data = readmatrix(obj.file_path, 'NumHeaderLines', Nheader);
               obj.data = data;
               save([obj.file_path(1:end-8), '.mat'], 'data');
            end
            toc
            
            % --- Data Preprocessing ---
            % Expand obj.indices for additional data columns
            obj.ind.axisSumPI = ind_cntr + (1:3);
            obj.ind.sinarg = obj.ind.debug(1);
            
            % Convert microseconds to seconds for time vector
            obj.time = (obj.data(:,obj.ind.time) - obj.data(1,obj.ind.time)) * 1.0e-6;
            
            % Unscale highResolutionGain
            if obj.para.blackbox_high_resolution
                blackbox_high_resolution_scale = 10.0;
                ind_bb_high_res = [obj.ind.gyroADC, obj.ind.gyroUnfilt, obj.ind.rcCommand, obj.ind.setpoint(1:3)];
                obj.data(:, ind_bb_high_res) = 1.0 / blackbox_high_resolution_scale * obj.data(:, ind_bb_high_res);
            end
            
            % Unscale and remap sinarg
            sinargScale = 5.0e3;
            obj.data(:,obj.ind.sinarg) = 1.0 / sinargScale * obj.data(:,obj.ind.sinarg);
                        
            % Create an additional entry for the pi sum
            obj.data = [obj.data, obj.data(:,obj.ind.axisP) + obj.data(:,obj.ind.axisI)];
            
            % Create different sampling times
            Ts      = obj.para.looptime * 1.0e-6;             % Gyro loop
            obj.Ts_cntr = obj.para.pid_process_denom * Ts;        % Control loop
            obj.Ts_log  = obj.para.frameIntervalPDenom * obj.Ts_cntr; % Logging loop
            
                       
         end

         function obj = plot_Gyro_Signal(obj, do_insert_legends)
            % --- Roll ---
            ax(1) = subplot(3,1,1); hold(ax(1),'on');
            plot(ax(1), obj.time, obj.data(:, obj.ind.setpoint(1)));
            plot(ax(1), obj.time, obj.data(:, obj.ind.gyroUnfilt(1)));
            plot(ax(1), obj.time, obj.data(:, obj.ind.gyroADC(1)));
            grid(ax(1),'on');
            ylabel(ax(1),'Roll (deg/s)');
            title(ax(1),'Gyro Signals - Setpoints from RemoteController');
        
            if do_insert_legends
                legend(ax(1), ...
                    {'setpoint','gyro unfiltert','gyroADC'}, ...
                    'Location','northeastoutside');
            end
        
            % --- Pitch ---
            ax(2) = subplot(3,1,2); hold(ax(2),'on');
            plot(ax(2), obj.time, obj.data(:, obj.ind.setpoint(2)));
            plot(ax(2), obj.time, obj.data(:, obj.ind.gyroUnfilt(2)));
            plot(ax(2), obj.time, obj.data(:, obj.ind.gyroADC(2)));
            grid(ax(2),'on');
            ylabel(ax(2),'Pitch (deg/s)');
            legend(ax(2),'off');
        
            % --- Yaw ---
            ax(3) = subplot(3,1,3); hold(ax(3),'on');
            plot(ax(3), obj.time, obj.data(:, obj.ind.setpoint(3)));
            plot(ax(3), obj.time, obj.data(:, obj.ind.gyroUnfilt(3)));
            plot(ax(3), obj.time, obj.data(:, obj.ind.gyroADC(3)));
            grid(ax(3),'on');
            ylabel(ax(3),'Yaw (deg/s)');
            xlabel(ax(3),'Time (s)');
            legend(ax(3),'off');
        
            % --- Axis linking & styling ---
            xmax = max(obj.time);
            linkaxes(ax,'x');
            xlim(ax,[0 xmax]);
            set(findall(gcf,'Type','line'),'LineWidth',obj.linewidth);
         end
     end
end
