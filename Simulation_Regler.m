%% Simulation Frequenzgang Drohen

clc, clear variables
addpath lib/
s = tf('s');

%% Laden von Drohnenflüge als Referenz
flight_folder = '20250908';

<<<<<<< HEAD:Simulation_Regler.m
% quad = 'aosmini';
% log_name = '20250907_aosmini_00.bbl.csv';
=======
quad = 'aosmini';
log_name = '20250908_aosmini_00.01.csv';
>>>>>>> efbc4e3b8be1b3f581d270aa6ac7f1a1890712ef:Simulation.m

quad = 'apex5';
log_name = '20250907_apex5_00.bbl.csv';

% quad = 'flipmini';
% log_name = '20250907_flipmini_00.bbl.csv';

% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;

%% Daten Laden, 

file_path = fullfile(flight_folder, log_name);
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Read the data
%  - If its the first time from the .csv and save a mat, otherwise the
%    .mat. This increases load speed significantly.
tic
try
   load([file_path(1:end-8), '.mat'])
catch exception
   data = readmatrix(file_path, 'NumHeaderLines', Nheader);
   save([file_path(1:end-8), '.mat'], 'data');
end
[Ndata, Nsig] = size(data);
toc
% Expand index
ind.axisSumPI = ind_cntr + (1:3);
ind.sinarg = ind.debug(1);

%% Einstellungen Bodeplot

% Defines
set(cstprefs.tbxprefs, 'MagnitudeUnits', 'dB');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'UnwrapPhase', 'Off');
set(cstprefs.tbxprefs, 'Grid', 'On');
set(groot,'defaultLineLineWidth',1.2);   % global für alle neuen Linien
set(0, 'defaultAxesColorOrder', get_my_colors);

% Bodeoptions
opt = bodeoptions('cstprefs');
opt.Xlim = { [0.1 1e3] };      % x-Achse: 0.1 Hz bis 1000 Hz

             

%% Filterparameter

%Sampletime
Ts = para.looptime * 1.0e-6;             % Gyro loop
z = tf('z', Ts);

%Drohnenparameter
switch quad
    case 'aosmini'
        % type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
        para_new.gyro_lpf            = 0;       % dono what this is
        para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
        para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
        para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
        para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
        para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
        para_new.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
        para_new.gyro_notch_cutoff   = [0, 0]; % 3 dB frequency where filter attenuates
        para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
        para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
        para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
        para_new.dterm_lpf2_hz       = 120;     % frequency of dterm lpf 2
        para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
        para_new.dterm_notch_hz      = 0;     % frequency of dterm notch
        para_new.dterm_notch_cutoff  = 0;     % 3 dB frequency where filter attenuates
        para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
        switch ind_ax
            case 1 % roll: [33, 52, 26, 0]
                P_new       = 33;
                I_ratio_new = 52/52;
                D_new       = 26;
            case 2 % pitch: [58, 98, 44, 0]
                P_new       = 58;
                I_ratio_new = 98/98;
                D_new       = 44;
            case 3 % yaw: [42, 65, 3, 0]
                P_new       = 42;
                I_ratio_new = 65/65;
                D_new       = 3;
        end
    case 'apex5'
        % type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
        para_new.gyro_lpf            = 0;       % dono what this is
        para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
        para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
        para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
        para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
        para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
        para_new.gyro_notch_hz       = [0, 520]; % frequency of gyro notch 1 and 2
        para_new.gyro_notch_cutoff   = [0, 448]; % 3 dB frequency where filter attenuates
        para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
        para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
        para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
        para_new.dterm_lpf2_hz       = 130;     % frequency of dterm lpf 2
        para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
        para_new.dterm_notch_hz      = 235;     % frequency of dterm notch
        para_new.dterm_notch_cutoff  = 202;     % 3 dB frequency where filter attenuates
        para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
        switch ind_ax
            case 1 % roll: [49, 83, 33, 0]
                P_new       = 49;
                I_ratio_new = 83/83;
                D_new       = 33;
            case 2 % pitch: [61, 103, 39, 0]
                P_new       = 61;
                I_ratio_new = 103/103;
                D_new       = 39;
            case 3 % yaw: [42, 104, 3, 0]
                P_new       = 42;
                I_ratio_new = 104/104;
                D_new       = 3;
        end
    case 'flipmini'
        % type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
        para_new.gyro_lpf            = 0;       % dono what this is
        para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
        para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
        para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
        para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
        para_new.gyro_soft2_type     = 0;       % type of gyro lpf 
        para_new.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
        para_new.gyro_notch_cutoff   = [0, 0]; % 3 dB frequency where filter attenuates
        para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
        para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
        para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
        para_new.dterm_lpf2_hz       = 140;     % frequency of dterm lpf 2
        para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
        para_new.dterm_notch_hz      = 0;       % frequency of dterm notch
        para_new.dterm_notch_cutoff  = 0;       % 3 dB frequency where filter attenuates
        para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
        switch ind_ax
            case 1 % roll: [46, 74, 30, 0]
                P_new       = 46;
                I_ratio_new = 74/74;
                D_new       = 30;
            case 2 % pitch: [71, 118, 47, 0]
                P_new       = 71;
                I_ratio_new = 118/118;
                D_new       = 47;
            case 3 % yaw: [35, 70, 3, 0]
                P_new       = 35;
                I_ratio_new = 70/70;
                D_new       = 3;
        end
    otherwise
        warning(' no valid quad selected');
end

%% Einlesen von alten PID Werten

% PID parameters
    pid_axis = {'rollPID', 'pitchPID', 'yawPID'};
    if (length(para.(pid_axis{ind_ax})) == 5)
        if (para.(pid_axis{ind_ax})(3) ~= para.(pid_axis{ind_ax})(4) && ...
                para.(pid_axis{ind_ax})(4) ~= 0)
            warning([pid_axis{ind_ax}, ' different D gains']);
        end
        % Remove dynamic D-Term
        para.(pid_axis{ind_ax}) = para.(pid_axis{ind_ax})([1 2 3 5]);
    end
    if para.(pid_axis{ind_ax})(4) ~= 0
        warning([pid_axis{ind_ax}, ' FF is not zero']);
    end
    % Insert 0 for FF
    PID = para.(pid_axis{ind_ax}) .* [get_pid_scale(ind_ax), 0];

%% PID Vektor

pid_scale = [get_pid_scale(ind_ax), 1];

PID_new(1) = P_new * pid_scale(1);                  % Proportionalanteil

%KI von KP abhängig, daher rückrechnung nötig
fI = PID(2) / (2 * pi * PID(1));        % alte Integralfrequenz aus altem PID
fI_new = fI * I_ratio_new;              % gewünschte skalierte Integralfrequenz

PID_new(2) = 2 * pi * PID_new(1) * fI_new;          % Integralanteil
%Differentialanteil unabhängig
PID_new(3) = D_new * pid_scale(3);                  % Differentialanteil
PID_new(4) = 0;

%% Funktion für Tiefpassfilter
function G = Tiefpassfilter(type, omega, s)
    if omega == 0
        G = tf(1);
    else
        switch type
            case 0 % PT1
                G = 1 / (1 + s/omega);
            case 1 % BIQUAD 
                G = 1 / (1 + s/omega); %unklar was für eine Übertrafungsfunktion BIQUAD ist
            case 2 % PT2
                G = 1 / (1 + s/omega)^2;
            case 3 % PT3
                G = 1 / (1 + s/omega)^3;
            otherwise
                error('Unbekannter Filtertyp');
        end
    end
end

%% Gyro-Tiefpassfilter berechnung

%Berechnung LPF1
omega1=2*pi*para_new.gyro_lowpass_hz;
G_LPF1 = Tiefpassfilter(para_new.gyro_soft_type, omega1, s);

%Berechnung LPF2
omega2 = 2*pi*para_new.gyro_lowpass2_hz;
G_LPF2 = Tiefpassfilter(para_new.gyro_soft2_type, omega2, s);

%Berechnung D-spezifischer Tiefpassfilter 1
omegad1 = 2*pi*para_new.dterm_lpf_hz;
G_LPFD1 = Tiefpassfilter(para_new.dterm_filter_type, omegad1, s);

%Berechnung D-spezifischer Tiefpassfilter 2
omegad2 = 2*pi*para_new.dterm_lpf2_hz;
G_LPFD2 = Tiefpassfilter(para_new.dterm_filter2_type, omegad2, s);

%% Funktionen diskretisierter Tiefpass
function G = lpf_pt1_discrete(fc, Ts)
    Om = 2*pi*fc*Ts;            % cut off frequency
    k  = Om/(1+Om);             % Steplenght
    G  = tf(k, [1 -(1-k)], Ts); %Transferfunction
end

function G = lpf_discrete(type, fc, Ts)
    if fc == 0      %in case LPF is not activatet
        G = tf(1,1,Ts);
        return;
    end
    switch type
        case 0  % PT1
            G = lpf_pt1_discrete(fc, Ts);

        case 2  % PT2 = PT1 ⨯ PT1, mit Cutoff-Correction
            c = 1.553773974;
            H1 = lpf_pt1_discrete(fc*c, Ts);
            G  = H1*H1;

        case 3  % PT3 = PT1 ⨯ PT1 ⨯ PT1, mit Cutoff-Correction
            c = 1.961459177;
            H1 = lpf_pt1_discrete(fc*c, Ts);
            G  = H1*H1*H1;

        case 1  % BIQUAD: ohne Q keine eindeutige Spec; Placeholder = PT1
            warning('BIQUAD ohne Q: ersetze provisorisch durch PT1.');
            G = lpf_pt1_discrete(fc, Ts);

        otherwise
            error('Unbekannter Filtertyp');
    end
end


%% Gyro-Lowpassfilter diskretisiert
%Calculation LPF1
G_LPF1_dis  = lpf_discrete(para_new.gyro_soft_type,  para_new.gyro_lowpass_hz,  Ts);
%Calculation LPF2
G_LPF2_dis  = lpf_discrete(para_new.gyro_soft2_type, para_new.gyro_lowpass2_hz, Ts);
%Calculation LPF1 D-share
G_LPFD1_dis = lpf_discrete(para_new.dterm_filter_type,  para_new.dterm_lpf_hz,  Ts);
%Calculation LPF2 D-share
G_LPFD2_dis = lpf_discrete(para_new.dterm_filter2_type, para_new.dterm_lpf2_hz, Ts);

%% Funktion für Notch-Filter
function G_notch = Notch(f0, fcut, s)
    if f0 ~= 0         
        f2 = f0^2 / fcut;
        Q  = f0 / (f2-fcut);           % Calculation with cutoff frequency
        omega = 2*pi*f0;
        G_notch = (s^2 + omega^2) / (s^2 + (omega/Q)*s + omega^2);     
    else
        G_notch = tf(1);
    end
end

%% Discretization of Notch filter
function G_notch_dis = Notch_dis(f0, fcut, Ts, z)
if f0 ~= 0         
        f2 = f0^2 / fcut;
        Q  = f0 / (f2-fcut);           % Calculation with cutoff frequency
        omega = 2*pi*f0*Ts;
        alpha = 1 - sin(omega)/2*Q;
       G_notch_dis = (1 - 2*cos(omega)*z^-1 + z^-2) ...
              / ((1 + alpha) - 2*cos(omega)*z^-1 + (1 - alpha)*z^-2);    
    else
        G_notch_dis = tf(1);
    end
end


%% Notch Filter berechnung

%Berechnung Notchfilter 1
G_Notch1 = Notch(para_new.gyro_notch_hz(1), para_new.gyro_notch_cutoff(1), s);
G_Notch1_dis = Notch_dis(para_new.gyro_notch_hz(1), para_new.gyro_notch_cutoff(1), Ts, z);

%Berechnung Notchfilter 2
G_Notch2 = Notch(para_new.gyro_notch_hz(2), para_new.gyro_notch_cutoff(2), s);
G_Notch2_dis = Notch_dis(para_new.gyro_notch_hz(2), para_new.gyro_notch_cutoff(2), Ts, z);

%Berechnung Notchfilter D-Anteil
G_NotchD = Notch(para_new.dterm_notch_hz, para_new.dterm_notch_cutoff, s);
G_NotchD_dis = Notch_dis(para_new.dterm_notch_hz, para_new.dterm_notch_cutoff, Ts, z);


%% Ausgabe PI und D Regler

C_PI = PID_new(1) + PID_new(2)/s;   %PI Regler ohne Filter
Kp = tf(PID_new(1));                %Für Simulink
C_PI_LFP_Notch = C_PI * G_LPF1 * G_LPF2 * G_Notch1 * G_Notch2;     %PI Regler mit Filter

C_D = PID_new(3)*s;     %D Regler ohne Filter
C_D_LFP = C_D*G_LPFD1*G_LPFD2;
C_D_LFP_Notch = C_D* G_LPFD1*G_LPFD2*G_NotchD;      %D Regler mit Filter

bereich = {2*pi*0.2, 2*pi*1000};

figure(1)
bode(C_PI_LFP_Notch,C_D_LFP_Notch, opt);
title('PI und D Regler');
legend('C_PI Sim','C_D Sim');


%% PID Vektor Diskretisiert

% Calculation of discretization of PI Part of the Controller

C_I_dis = (PID_new(2)*Ts) / (1 - z^-1);   % Calculation I Part 
C_PI_dis = PID_new(1) + C_I_dis;          % Combination P and I Part of Controller
C_PI_LFP_Notch_dis = C_PI_dis * G_LPF1_dis * G_LPF2_dis;     %PI Controller with LPF

%Calculation of D Part of the Controller
C_D_dis = (PID_new(3)*(1 - z^-1))/Ts;     %Calculation D Part
C_D_LFP_dis = C_D_dis*G_LPFD1_dis*G_LPFD2_dis;  %D Part with LPF
figure(7)
bode(C_PI_LFP_Notch_dis, C_D_LFP_dis, C_PI_LFP_Notch, C_D_LFP_Notch, opt);
title('PI und D Regler');


%% Simulation Strecke P

switch quad
    case 'aosmini' 
        switch ind_ax
            case 1  %Roll
                %Relevante Frequenzen
                 w2 = 11*2*pi();
                 w1 = 14*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
            case 2  %Pitch
                %Relevante Frequenzen
                 w2 = 10*2*pi();
                 w1 = 13*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
            case 3  %Yaw
                %Relevante Frequenzen
                 w2 = 60*2*pi();
                 w1 = 10*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
        end

    case 'apex5' 
        switch ind_ax
            case 1  %Roll
                %Relevante Frequenzen
                 w2 = 10*2*pi();
                 w1 = 13*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion

            case 2  %Pitch
                %Relevante Frequenzen
                 w2 = 10*2*pi();
                 w1 = 13*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
            case 3  %Yaw
                %Relevante Frequenzen
                 w2 = 60*2*pi();
                 w1 = 10*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
        end

        case 'flipmini' 
        switch ind_ax
            case 1  %Roll
                %Relevante Frequenzen
                 w2 = 10*2*pi();
                 w1 = 12*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion

            case 2  %Pitch
                %Relevante Frequenzen
                 w2 = 10*2*pi();
                 w1 = 13*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
            case 3  %Yaw
                %Relevante Frequenzen
                 w2 = 60*2*pi();
                 w1 = 10*2*pi();
                 wt = 60*2*pi();
                        
                 %Übertragungsfunktionen abgeschätzt
                 G1 = 1 / (s/w1);
                 G2 = 1 / (1 + (s/w2));
                 Gt = exp(-s * (1/wt));     %Totzeit
                 P_ges = G1*G2*Gt;          %Gesamtübertragungsfunktion
        end
end

%% Plot P
figure(3)
bode(P_ges);
title('Übertragungsfunktion P');

%% Step response Drone

