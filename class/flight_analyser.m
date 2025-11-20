classdef flight_analyser
    properties       
        % Raw flight data inputs

        time             % Time vector [s]
        data             % Flgiht data
        ind              % Columns
        Ts_log           % Logging Time [s]

        % Used in different functions
        freq_spectra
        spectra
        freq_spectogram
        spectrogram_unf
        spectrogram_fil
        throttle_all

        % Default values
        linewidth            (1,1) double  = 1.2
        
        
    end

    methods
        function obj = flight_analyser(data, ind, Ts_log)
            obj.data = data;
            obj.ind = ind;
            obj.Ts_log = Ts_log;

        end

        function obj = calculate_spectra(obj, Nestfaspec, koverlapspec)
            % =============================================================
            %  Spectral Analysis
            % =============================================================
            % Performs frequency domain analysis of gyro signals and PID sums
            % 
            % Steps:
            % 1. Configure spectral estimation parameters
            % 2. Calculate power spectra using Hann window
            % 3. Convert power to amplitude spectra
            % 4. Store results for unfiltered, filtered gyro and axis sums
                
            
            % Select data columns for spectral analysis
            data_for_spectra = obj.data(:,[obj.ind.gyroUnfilt, ...
                                       obj.ind.gyroADC, ...
                                       obj.ind.axisSum, ...
                                       obj.ind.setpoint(1:3)]);
            
            % Configure spectral estimation parameters
            Nestspec = round(Nestfaspec / obj.Ts_log);    % Window length in samples
            windowspec   = hann(Nestspec, 'periodic');        % Hann window for better frequency resolution
            Noverlap = floor(koverlapspec * Nestspec);    % Overlap for smoother estimates

            % Calculate power spectral density
            [pxx, obj.freq_spectra] = estimate_spectra(data_for_spectra, windowspec, Noverlap, Nestspec, obj.Ts_log);
            obj.spectra = sqrt(pxx); % Convert power to amplitude spectra (dc needs to be scaled differently)


        end
        function obj = calculate_spectogram(obj, Nestfaspec, koverlapspec)

            Nestspec = round(Nestfaspec / obj.Ts_log);    % Window length in samples
            windowspec   = hann(Nestspec, 'periodic');        % Hann window for better frequency resolution
            Noverlap = floor(koverlapspec * Nestspec);    % Overlap for smoother estimates

            % Calculate throttle resolution
            Nres = floor(max(obj.data(:,obj.ind.setpoint(4))) / 1e1 / 2); % should give 40 at 80% throttle constrain

            % Initialize storage for multiple spectrograms
            num_spectrograms = 3;    % One per axis
            obj.spectrogram_unf = cell(1, num_spectrograms);
            obj.spectrogram_fil = cell(1, num_spectrograms);
            obj.freq_spectogram = cell(1, num_spectrograms);
            obj.throttle_all = cell(1, num_spectrograms);

            % Calculate spectrograms for unfiltered gyro data
            for spectrogram_nr = 1:num_spectrograms
                [pxx, freq, throttle] = estimate_spectrogram( ...
                    obj.data(:, obj.ind.gyroUnfilt(spectrogram_nr)), ...    % Raw gyro data
                    obj.data(:, obj.ind.setpoint(4)) / 10.0, ...            % Throttle values
                    windowspec, Noverlap, Nestspec, Nres, obj.Ts_log);

                obj.spectrogram_unf{spectrogram_nr} = sqrt(pxx); % Convert power to amplitude
                obj.freq_spectogram{spectrogram_nr} = freq;
                obj.throttle_all{spectrogram_nr} = throttle;
            end

                 % Calculate spectrograms for filtered gyro data
            for spectrogram_nr = 1:num_spectrograms
                [pxx, freq, throttle] = estimate_spectrogram( ...
                    obj.data(:, obj.ind.gyroADC(spectrogram_nr)), ...    % Filtered gyro data
                    obj.data(:, obj.ind.setpoint(4)) / 10.0, ...         % Throttle values
                    windowspec, Noverlap, Nestspec, Nres, obj.Ts_log);

                obj.spectrogram_fil{spectrogram_nr} = sqrt(pxx); % Convert power to amplitude
                obj.freq_spectogram{spectrogram_nr} = freq;
                obj.throttle_all{spectrogram_nr} = throttle;
            end
            
        end

        function plot_spectra(obj, do_insert_legends)
             % =================================================================
            %  FIGURE GYRO SPECTRA
            % =================================================================
            % Generates spectral analysis plots
      
            ax(1) = subplot(3, 1, 1);
            plot(ax(1), obj.freq_spectra, obj.spectra(:,1:3));
            grid(ax(1),'on');
            ylabel(ax(1),'Gyro (deg/s)');
            set(ax(1),'YScale','log');
            title(ax(1),'Unfiltered gyro magnitude spectra');
            legend(ax(1), 'off');
            
            % ---------- Subplot 2: Filtered (ADC) gyro spectra ----------
            ax(2) = subplot(3, 1, 2);
            plot(ax(2), obj.freq_spectra, obj.spectra(:,4:6));
            grid(ax(2),'on');
            ylabel(ax(2),'Gyro (deg/s)');
            set(ax(2),'YScale','log');
            title(ax(2),'Filtered (ADC) gyro magnitude spectra');
            
            if do_insert_legends
                legend(ax(2), {'Roll','Pitch','Yaw'}, 'Location','northeast');
            end
            
            % ---------- Subplot 3: Axis sum spectra ----------
            ax(3) = subplot(3, 1, 3);
            plot(ax(3), obj.freq_spectra, obj.spectra(:,7:9));
            grid(ax(3),'on');
            ylabel(ax(3),'AxisSum');
            xlabel(ax(3),'Frequency (Hz)');
            set(ax(3),'YScale','log');
            title(ax(3),'Axis sum spectra');
            legend(ax(3), 'off');
            
            % ---------- Link axes and set Nyquist limit ----------
            linkaxes(ax, 'x');
            
            nyq = 1/(2 * obj.Ts_log);   % Nyquist frequency
            xlim([0, nyq]);
            
            % Optional uniform y-limits
            try
                ylim(ax, [1e-3 1e1]);
            catch
            end
            
            % Global line width
            set(findall(gcf,'type','line'), 'LineWidth', obj.linewidth);

        end
        function plot_spectogram(obj, num_spectrograms)
              % =================================================================
            %  FIGURE GYRO SPECTRA
            % =================================================================
            % Generates spectogram plots
                 

            sgtitle('Gyro Spectrograms')
            axes_labels = {'Roll', 'Pitch', 'Yaw'};
            c_lim = [5e-2 3e0];
            colormap('jet')

           % --- Unfiltered Spectrograms ---
            for spectrogram_nr = 1:num_spectrograms
                subplot(230 + spectrogram_nr)
                qmesh = pcolor(obj.freq_spectogram{spectrogram_nr}, ...
                               obj.throttle_all{spectrogram_nr}, ...
                               obj.spectrogram_unf{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
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
                qmesh = pcolor(obj.freq_spectogram{spectrogram_nr}, ...
                               obj.throttle_all{spectrogram_nr}, ...
                               obj.spectrogram_fil{spectrogram_nr});
                set(qmesh, 'EdgeColor', 'none');
                xlabel('Frequency (Hz)')
                if spectrogram_nr == 1
                    ylabel('Throttle (%)')
                end
                title([axes_labels{spectrogram_nr}, ' – with Filter'])
                colormap('jet')
                set(gca, 'ColorScale', 'log')
                clim(c_lim);
                ylim([0 100])
            end
        end
    end

end




