function splot(mystring)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Alternative zu subplot(211) und subplot(212)
% M.E.Peter, 14.12.2011
%
% splot(strin)
% -----------
%
% string: '21'  -> 1ster Plot Zeitbereich
% string: '22'  -> 2ter  Plot Zeitbereich
% string: '21b'  -> 1ster Plot Frequenzbereich
% string: '22b'  -> 2ter  Plot Frequenzbereich
%
% Bsp:
% figure(1)
% splot('11b')
% semilogx(w,db(abs(ak(:)))),grid on
% ylabel('Amplitude [dB]')
% splot('12b')
% semilogx(w,p(:)),grid on
% ylabel('Phase [°]')
% xlabel('Frequenz [rad/s]')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mystring
        % 50/50 Aufteilung
        case '21'
            vec = [0.1300    0.5838-0.025    0.7750    0.3412+0.025];           
        case '22'
            vec = [0.1300    0.1100    0.7750    0.3412+0.025];
        % 60/40 Aufteilung
        case '21b'
            vec = [0.1300    0.6097-0.15    0.7750    0.3153+0.15];           
        case '22b'
            set(gca,'xticklabel',{})
            vec = [0.1300    0.1358    0.7750    0.3153-0.05];
        case '21bt'
            vec = [0.1300    0.6097-0.15    0.7750    0.3153+0.15];           
        case '22bt'
            vec = [0.1300    0.1358    0.7750    0.3153-0.05];
        otherwise
                return;
end
subplot('position',vec); 

return
