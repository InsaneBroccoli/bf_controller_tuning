classdef ControllerAnalyzer
    % Berechnet alte und neue Closed-Loop-Systeme (T, S, SC, SP)

    properties
        freq_resp_result
        Cpi_ana_new
        Cd_ana_new
        Gf_ana_new
        Ts_log
        do_compensate_iterm
    end

    methods
        function obj = ControllerAnalyzer(freq_resp_result, Cpi_ana_new, Cd_ana_new, Gf_ana_new, Ts_log, do_compensate_iterm)
            obj.freq_resp_result = freq_resp_result;
            obj.Cpi_ana_new = Cpi_ana_new;
            obj.Cd_ana_new = Cd_ana_new;
            obj.Gf_ana_new = Gf_ana_new;
            obj.Ts_log = Ts_log;
            obj.do_compensate_iterm = do_compensate_iterm;
        end

        function [CL_ana, CL_ana_new] = calculate_closed_loops(obj)
            % Alte und neue Closed-Loop-Systeme berechnen
            CL_ana = controller_analysis.calculate_closed_loop( ...
                obj.freq_resp_result.Cpi_ana, tf(1,1,obj.Ts_log), ...
                obj.freq_resp_result.P / obj.freq_resp_result.Gf_ana, ...
                obj.freq_resp_result.Gf_ana, obj.freq_resp_result.Cd_ana);

            CL_ana_new = controller_analysis.calculate_closed_loop( ...
                obj.Cpi_ana_new, tf(1,1,obj.Ts_log), ...
                obj.freq_resp_result.P / obj.freq_resp_result.Gf_ana, ...
                obj.Gf_ana_new, obj.Cd_ana_new);

            % Optionale Kompensation des I-Teils
            if obj.do_compensate_iterm
                Cpi_com = obj.freq_resp_result.Cpi / obj.freq_resp_result.Cpi_ana;
                CL_ana_ = calculate_closed_loop( ...
                    obj.freq_resp_result.Cpi_ana * Cpi_com, ...
                    tf(1,1,obj.Ts_log), ...
                    obj.freq_resp_result.P / obj.freq_resp_result.Gf_ana, ...
                    obj.freq_resp_result.Gf_ana, obj.freq_resp_result.Cd_ana);
                CL_ana_new_ = calculate_closed_loop( ...
                    obj.Cpi_ana_new * Cpi_com, ...
                    tf(1,1,obj.Ts_log), ...
                    obj.freq_resp_result.P / obj.freq_resp_result.Gf_ana, ...
                    obj.Gf_ana_new, obj.Cd_ana_new);
                CL_ana.T = CL_ana_.T;
                CL_ana_new.T = CL_ana_new_.T;
            end
        end
    end
end
