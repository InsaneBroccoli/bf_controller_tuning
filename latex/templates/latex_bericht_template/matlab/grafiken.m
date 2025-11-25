%% Exportieren von Grafiken mit pdfexport.m
% pmic, 14.10.2013

clc,clear all
%%

x = rand(1e2,1);

figure(1)
stairs(x),grid on
xlabel('k')
ylabel('y = f(x)')
legend('y',4)
pdfexport('grafik1',[9.5 8.5])

figure(2)
splot('21')
stairs(x),grid on
ylabel('y = f(x)')
legend('y',4)
splot('22')
stairs(x),grid on
xlabel('k')
ylabel('y = f(x)')
legend('y',4)
pdfexport('grafik2') % -> [15 10]

G = tf(1,[1 1]);
f = logspace(-3,2,1e3);
g = squeeze(freqresp(G,2*pi*f));

figure(3)
splot('21b')
semilogx(f,db(abs(g))),grid on,hold on
xlim([f(1) f(end)])
ylabel('Amplitude [dB]')
splot('22b')
semilogx(f,angle(g)*180/pi),grid on,hold on
xlabel('Frequenz [Hz]')
xlim([f(1) f(end)])
ylabel('Phase [°]')
set(gca,'ytick',-180:90:180)
pdfexport('grafik3',[14 9])