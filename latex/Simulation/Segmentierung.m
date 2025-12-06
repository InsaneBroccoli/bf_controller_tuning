% Beispiel: Durchgehendes Signal in Segmente aufteilen

% --- Signal erzeugen ---
fs = 1000;                   % Samplingfrequenz [Hz]
t = 0 : 1/fs : 1;            % 1 Sekunde
x = sin(2*pi*50*t) + 0.2*randn(size(t));   % Sinus + Rauschen

% --- Parameter ---
segmentLen = 200;           % Anzahl Samples pro Segment

% --- Segmentierung ---
numSegments = floor(length(x) / segmentLen);

segments = cell(numSegments,1);
for k = 1:numSegments
    idx = (k-1)*segmentLen + 1 : k*segmentLen;
    segments{k} = x(idx);
end

% --- Beispiel: erstes Segment plotten ---

figure(1)
% --- Segment 1 ---
idx1 = 1 : segmentLen;
subplot(3,1,1)
plot(t(idx1), segments{1});
title('Segment 1');
xlabel('Zeit [s]');
ylabel('Amplitude');

% --- Segment 2 ---
idx2 = segmentLen+1 : 2*segmentLen;
subplot(3,1,2)
plot(t(idx2), segments{2});
title('Segment 2');
xlabel('Zeit [s]');
ylabel('Amplitude');

% --- Segment 3 ---
idx3 = 2*segmentLen+1 : 3*segmentLen;
subplot(3,1,3)
plot(t(idx3), segments{3});
title('Segment 3');
xlabel('Zeit [s]');
ylabel('Amplitude');
