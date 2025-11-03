classdef main_class
    properties
        file_path
        para_new
        ind_ax
        do_compensate_iterm 
        P_new
        I_ratio_new
        D_new
        unfgyroData
        gyroData
        setpoint
        motorData
        time
        unfgyroSpec
        adcgyroSpec
        axisSumSpec
        Specfreq
        Ts_log
        axisSumData
        transfData
        transfCoher
        transfOmega
        transfCpi
        transfCD
        transfCpiAna
        transfCDAna
        CloLoAan
        CloLoAanNew
        transfT
        step_resp_tra
        step_resp_com
        step_resp_tim
        SpecAmplUnf
        SpecfreqUnf
        SpecthroUnf
        SpecAmplFil
        SpecthroFil
        SpecfreqFil
        num_spectrograms
        Nestfaspec
        koverlapspec
        Nestfatra
        koverlaptra
    end

    methods
        % =================================================================
        %  ???
        % =================================================================
        function obj = main_class(file_path, para_new, ind_ax, do_compensate_iterm, ...
                P_new, I_ratio_new, D_new, Nestfaspec,koverlapspec,Nestfatra, koverlaptra)
            % Save all inputs
            obj.file_path = file_path;
            obj.para_new = para_new;
            obj.ind_ax = ind_ax;
            obj.do_compensate_iterm = do_compensate_iterm;
            obj.P_new = P_new;
            obj.I_ratio_new = I_ratio_new;
            obj.D_new = D_new;
            obj.Nestfaspec = Nestfaspec;
            obj.koverlapspec = koverlapspec;
            obj.Nestfatra = Nestfatra;
            obj.koverlaptra = koverlaptra;
        end

        % =================================================================
        %  ???
        % =================================================================
        function obj = run(obj)
           % Extract header information
            [para, Nheader, ind, ind_cntr] = extract_header_information(obj.file_path);
            
            % Read the data
            %  - If its the first time from the .csv and save a mat, otherwise the
            %    .mat. This increases load speed significantly.
            tic
            try
               load([obj.file_path(1:end-8), '.mat'])
            catch exception
               data = readmatrix(obj.file_path, 'NumHeaderLines', Nheader);
               save([obj.file_path(1:end-8), '.mat'], 'data');
            end
            toc
            
            % Expand index
            ind.axisSumPI = ind_cntr + (1:3);
            ind.sinarg = ind.debug(1);
            
            % Convert and evaluate time
            obj.time = (data(:,ind.time) - data(1,ind.time)) * 1.0e-6;
            
            % Unscale highResolutionGain
            if para.blackbox_high_resolution
                blackbox_high_resolution_scale = 10.0;
                ind_bb_high_res = [ind.gyroADC, ind.gyroUnfilt, ind.rcCommand, ind.setpoint(1:3)];
                data(:, ind_bb_high_res) = 1.0 / blackbox_high_resolution_scale * data(:, ind_bb_high_res);
            end
            
            % Unscale and remap sinarg
            sinargScale = 5.0e3;
            data(:,ind.sinarg) = 1.0 / sinargScale * data(:,ind.sinarg);
            
            % Assign negative sign for pid error
            data(:,ind.axisError) = -data(:,ind.axisError);
            
            % Create an additional entry for the pi sum
            data = [data, data(:,ind.axisP) + data(:,ind.axisI)];
            
            % Create different sampling times
            Ts      = para.looptime * 1.0e-6;             % Gyro loop
            Ts_cntr = para.pid_process_denom * Ts;        % Control loop
            obj.Ts_log  = para.frameIntervalPDenom * Ts_cntr; % Logging loop
            
            % Get evaluation index where Chirp was active
            ind_eval = get_ind_eval(data(:,ind.sinarg), data(:,ind.gyroADC(obj.ind_ax)));
            data(~ind_eval,ind.sinarg) = 0.0;
            T_eval_tot = size(data(ind_eval,ind.sinarg), 1) * obj.Ts_log
            
            % Calculate average throttle
            throttle_avg = median(data(ind_eval,ind.setpoint(4))) / 1.0e3;
            
            
            %% 
            % =============================================================
            %  show Gyro to select Teval and spectra (gyro and pid sum)
            % =============================================================
                
            obj.setpoint    = data(:, ind.setpoint(1:4));
            obj.unfgyroData = data(:, ind.gyroUnfilt(1:3));
            obj.gyroData    = data(:, ind.gyroADC(1:3));
            obj.motorData = data(:,ind.motor);
            obj.axisSumData = data(:, ind.axisSum);
                    
            % Select data for spectra
            data_for_spectra = data(:,[ind.gyroUnfilt, ...
                                       ind.gyroADC, ...
                                       ind.axisSum, ...
                                       ind.setpoint(1:3)]);
            
            % Parameters
            Nestspec = round(obj.Nestfaspec / obj.Ts_log);
            windowspec   = hann(Nestspec, 'periodic');
            Noverlap = floor(obj.koverlapspec * Nestspec);
            [pxx, freq] = estimate_spectra(data_for_spectra, windowspec, Noverlap, Nestspec, obj.Ts_log);
            spectra = sqrt(pxx); % power -> amplitude (dc needs to be scaled differently)

            % Output Frequencystep
            obj.Specfreq = freq;
            
            % Output unfiltert Gyro Spectra
            obj.unfgyroSpec = spectra(:,1:3);
       
            % Output filtert Gyro Spectra
            obj.adcgyroSpec = spectra(:,4:6);
            
            % Output AxisSum Spectra
            obj.axisSumSpec = spectra(:,7:9);                     
            
        %%
        % =============================================================
        %  Spectrogram
        % =============================================================
            % Parameters

            Nres     = floor(max(data(:,ind.setpoint(4))) / 1e1 / 2) % should give 40 at 80% throttle constrain

        
            % Initialisierung
            obj.num_spectrograms = 3;
            spectrograms = cell(1, obj.num_spectrograms);
            freq_all = cell(1, obj.num_spectrograms);
            throttle_all = cell(1, obj.num_spectrograms);
            
            % --- Berechnung der Spektrogramme im Loop ---
            for spectrogram_nr = 1:obj.num_spectrograms
                [pxx, freq, throttle] = estimate_spectrogram( ...
                    data(:, ind.gyroUnfilt(spectrogram_nr)), ...
                    data(:, ind.setpoint(4)) / 10.0, ...
                    windowspec, Noverlap, Nestspec, Nres, obj.Ts_log);
            
                spectrograms{spectrogram_nr} = sqrt(pxx); % power -> amplitude
                freq_all{spectrogram_nr} = freq;
                throttle_all{spectrogram_nr} = throttle;
            end
            
            obj.SpecAmplUnf = spectrograms;
            obj.SpecfreqUnf = freq_all;
            obj.SpecthroUnf = throttle_all;
           
            for spectrogram_nr = 1:obj.num_spectrograms
                [pxx, freq, throttle] = estimate_spectrogram( ...
                    data(:, ind.gyroADC(spectrogram_nr)), ...
                    data(:, ind.setpoint(4)) / 10.0, ...
                    windowspec, Noverlap, Nestspec, Nres, obj.Ts_log);
            
                spectrograms{spectrogram_nr} = sqrt(pxx); % power -> amplitude
                freq_all{spectrogram_nr} = freq;
                throttle_all{spectrogram_nr} = throttle;
            end
            obj.SpecAmplFil = spectrograms;
            obj.SpecfreqFil = freq_all;
            obj.SpecthroFil = throttle_all; 
         
            %% 
            % =============================================================
            %  Frequency response estimation and calculation
            % =============================================================

            % Parameters
            Nest     = round(obj.Nestfatra / obj.Ts_log);
            Noverlap = floor(obj.koverlaptra * Nest);
            window   = hann(Nest, 'periodic');
            
            % Linear filter for zero phase excitation filter (apply_rotfiltfilt)
            Dlp = sqrt(3) / 2;
            wlp = 2 * pi * 10;
            Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), obj.Ts_log, 'tustin');
            
            % T  , Gyw: w -> y
            inp = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.setpoint(obj.ind_ax)));
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.gyroADC(obj.ind_ax)) );
            [T, C_T] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, obj.Ts_log);
            
            % SCw, Guw: w -> u
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSum(obj.ind_ax)));
            [Guw, C_Guw] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, obj.Ts_log);
            
            %      Gvw: w -> v (v := u only from PI cntrl)
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSumPI(obj.ind_ax)));
            [Gvw, C_Gvw] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, obj.Ts_log);
            
            % P  , Gyu: u -> y
            P = T / Guw;
            
            % % P  , Gyu: u -> y (direct measurement, results are slightly worse)
            % inp = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSum(ind_ax)));
            % out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.gyroADC(ind_ax)));
            % [Pd, C_Pd] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, Ts_log);
            
            % Calculated controller frequency response estimates
            Cpi = Gvw / (1 - T);
            Cd  = Guw * Gvw / T * (1 / Guw - 1 / Gvw);
            
            % Index and frequency for bode plots
            omega_bode = 2*pi*P.Frequency;
            
            
            % Downsample analytical controller transferfunction and convert to frd objects
            
            [Cpi_ana, Cd_ana, Gf_ana, PID, para_used] = ...
                calculate_transfer_functions(para, obj.ind_ax, throttle_avg, Ts_cntr);
            
            if Gf_ana.Ts < obj.Ts_log % by using Gf_ana.Ts we secure that we do this only once
                Gf_ana  = downsample_frd(Gf_ana , obj.Ts_log, P.Frequency);
                Cpi_ana = downsample_frd(Cpi_ana, obj.Ts_log, P.Frequency);
                Cd_ana  = downsample_frd(Cd_ana , obj.Ts_log, P.Frequency);
            end
            
            % =============================================================
            %  ???
            % =============================================================
            
            
            %% Plant and used controllers
            obj.transfData = P / Gf_ana;
            obj.transfCoher = C_T * C_Guw;
            obj.transfOmega = omega_bode;
            obj.transfCpi = Cpi;
            obj.transfCD = Cd;
            obj.transfCpiAna = Cpi_ana;
            obj.transfCDAna = Cd_ana;
            obj.transfT = T;
            
            %% New controller and filter parameters
            
            tic
            
            pid_axis = {'rollPID', 'pitchPID', 'yawPID'};
            
            % PID parameters
            fprintf('   used PID parameters are:\n');
            fprintf(['      ', pid_axis{obj.ind_ax}, ':  %d, %d, %d\n'], ...
                para.(pid_axis{obj.ind_ax})(1:3));
            
            % Inform user about parameters
            para_used_fieldnames = fieldnames(para_used);
            Npara_used = size(para_used_fieldnames, 1);
            fprintf('   used parameters are:\n');
            for i = 1:Npara_used
                fprintf(['      ', para_used_fieldnames{i},': %d\n'], eval(['round(', 'para_used.', para_used_fieldnames{i}, ');']));
            end
            
            % Scale to new PID parameters
            pid_scale = [get_pid_scale(obj.ind_ax), 1];
            PID_new(1) = obj.P_new * pid_scale(1);
            fI         = PID(2) / (2 * pi * PID(1)); % extract fn from initial parametrization
            fI_new     = fI * obj.I_ratio_new;
            PID_new(2) = 2 * pi * PID_new(1) * fI_new;
            PID_new(3) = obj.D_new * pid_scale(3);
            PID_new(4) = 0;
            
            fprintf('   used fI is: %0.2f Hz\n\n', fI);
            
            % New PID parameters
            fprintf('   new PID parameters are:\n');
            obj.para_new.(pid_axis{obj.ind_ax}) = round( PID_new ./ pid_scale);
            obj.para_new.(pid_axis{obj.ind_ax}) = [obj.para_new.(pid_axis{obj.ind_ax})(1:3), ...
                                           obj.para_new.(pid_axis{obj.ind_ax})(3), ...
                                           obj.para_new.(pid_axis{obj.ind_ax})(4)];
            fprintf(['      ', pid_axis{obj.ind_ax}, ':  %d, %d, %d\n'], ...
                obj.para_new.(pid_axis{obj.ind_ax})(1:3));
            
            [Cpi_ana_new, Cd_ana_new, Gf_ana_new, PID_new, para_used_new] = ...
                calculate_transfer_functions(obj.para_new, obj.ind_ax, throttle_avg, Ts_cntr);
            
            % Inform user about new parameters
            para_used_fieldnames_new = fieldnames(para_used_new);
            Npara_used_new = size(para_used_fieldnames_new, 1);
            fprintf('   new parameters are:\n');
            for i = 1:Npara_used_new
                fprintf(['      ', para_used_fieldnames_new{i},': %d\n'], ...
                    eval(['round(', 'para_used_new.', para_used_fieldnames_new{i}, ');']));
            end
            
            fprintf('   new used fI is: %0.2f Hz\n\n', fI_new);
            
            % Downsample analytical controller transferfunction and convert to frd objects
            if Gf_ana_new.Ts < obj.Ts_log % by using Gf_ana.Ts we secure that we do this only once
                Gf_ana_new  = downsample_frd(Gf_ana_new , obj.Ts_log, P.Frequency);
                Cpi_ana_new = downsample_frd(Cpi_ana_new, obj.Ts_log, P.Frequency);
                Cd_ana_new  = downsample_frd(Cd_ana_new , obj.Ts_log, P.Frequency);
            end
            
            CL_ana     = calculate_closed_loop(Cpi_ana    , tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana    , Cd_ana    );
            CL_ana_new = calculate_closed_loop(Cpi_ana_new, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana_new, Cd_ana_new);
            if obj.do_compensate_iterm
                % Compensate only PI part
                Cpi_com = Cpi / Cpi_ana;
                CL_ana_      = calculate_closed_loop(Cpi_ana     * Cpi_com, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana    , Cd_ana    );
                CL_ana_new_  = calculate_closed_loop(Cpi_ana_new * Cpi_com, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana_new, Cd_ana_new);
                CL_ana.T     = CL_ana_.T;
                CL_ana_new.T = CL_ana_new_.T;
            end

            obj.CloLoAan = CL_ana;
            obj.CloLoAanNew = CL_ana_new;
            
  
            
            % Step responses
            f_max = min([para.dyn_notch_min_hz, para.gyro_rpm_notch_min]);
            T_mean = 0.1 * [-1, 1] + (Nest * obj.Ts_log) / 2;
            step_time = (0:Nest-1).'*obj.Ts_log;
            
            % Actual controller parameters
            step_resp = [calculate_step_response_from_frd(CL_ana.T    , f_max), ...
                         calculate_step_response_from_frd(CL_ana_new.T, f_max), ...
                         calculate_step_response_from_frd(T           , f_max)];
            step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
            step_resp = step_resp ./ step_resp_mean;

            obj.step_resp_tra = step_resp;
            obj.step_resp_tim = step_time;
            
            % New controller parameters
            step_resp = [calculate_step_response_from_frd(CL_ana.SP    , f_max), ...
                         calculate_step_response_from_frd(CL_ana_new.SP, f_max)];
            step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
            step_resp = step_resp - step_resp_mean;

            obj.step_resp_com = step_resp;
                       
            toc
        end
    end
end