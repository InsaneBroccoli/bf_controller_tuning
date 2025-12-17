# Projektarbeit

## Allgemein

### Titel

Betaflight - Offline Controller Tuning

### Beschreibung

Betaflight ist ein Open-Source Projekt, das sich auf die Regelung von FPV Quadcopter spezialisiert hat, dies insbesondere im Bereich von Race-Dronen. Die Firmware wird von einer grossen Community weiterentwickelt und ist auf vielen verschiedenen Hardware-Plattformen verfügbar. Betaflight ist die am häufigsten verwendete Firmware für FPV-Drohnen und wird von vielen Hobbyisten und Profis genutzt.

Für die Flugperformace ist die Qualität der Raten-Regelung von entscheidender Bedeutung. Die Regler sind dafür verantwortlich, dass die Drohne stabil fliegt, Störungen wie Wind und Strömungsabrisse (Prop-Wash) unterdrückt und möglichst agil auf Steuerbefehle reagiert. Die Regler und Filter müssen daher sorgfältig getunt und abgestimmt werden, um eine optimale Flugperformance zu erzielen. FPV Piloten fliegen üblicherweise im ACRO mode, sprich das System ist Ratengeregelt und die Stickinputs des Piloten werden als Sollwerte auf den Ratenregler gegeben. Das Tuning dieser Regler geschieht in der Community im Allgemeinen durch Trial and Error. Das bedeutet, dass die Piloten die Reglerparameter anpassen und dann testen. Dies kann jedoch zeitaufwendig und frustrierend sein, da es oft schwierig ist, insbesondere für Leien, die richtigen Parameter zu finden.

Bereits bestehend ist ein Proof of Concept, welches es ermöglicht die Raten-Regler offline basierend auf einem Messdaten-Modell der Drone zu tunen (basierend auf Chirp-Signal Experiment während dem Flug). Das Projekt soll auf diesem Proof of Concept aufbauen und die Methode des Offline Tunings für die Raten-Regler (Gyro) weiterentwickeln. Ziel ist es, eine vollständige, einfach erweiterbare und sauber implementierte Open-Source Libary zu entwickeln. Dies sowohl in MATLAB als auch in Python. Ausserdem soll für die Python Implementation eine ausführliche Dokumentation mit Beispielen für Anwender erstellt werden (Github Markdown Dokumentation). Das Ziel ist es, dass die Libary von der Community genutzt werden kann. Dies soll dazu beitragen, die Flugperformance der Drohnen zu verbessern und den Piloten eine einfachere Möglichkeit zu bieten, ihre Drohnen zu tunen.

### Studenten

Bianchi Yuri (biancyur) <br>
Dort Janick (dortjan1) <br>
Jurietti Dario (juriedar)

### Betreuer

Hauptbetreuer: Michael Peter (pmic) <br>
Nebenbetreuer: Prof. Dr. Ruprecht Altenburger (altb)

### Weitere Personen

Camille Huber (hurc)

### Voraussetzungen

- GRT, RT1
- Zumindest eine Person ist im Besitz des Fernpiloten-Zeugnis A1/A3 (A2 optional jedoch empfohlen)
- Genügend hohe Haftpflichtversicherung, muss von den Studierenden selbst organisiert werden
- Erfahrung mit FPV Dronen (Bauen, Reparieren, Konfigurieren, Fliegen)
- MATLAB, Python

### Vereinbarung

Die Vereinbarung wird zwischen dem Betreuer Michael Peter (pmic) und den Studierenden Yuri Bianchi (biancyur), Janick Dort (dortjan1) und Dario Jurietti (juriedar).

Da es sich um eine dreier PA handelt wird das Projekt in drei separat bewertbare Teile aufteilt:
- MATLAB Implementation
- Python Implementation
- Markdown Dokumentation mit Beispielen der Matlab und Python Implementation

### Zeitplan

Siehe: https://intra.zhaw.ch/departemente/school-of-engineering/bachelorstudium/projekt-und-bachelorarbeiten


| Projektarbeit HS 25          |                                                              |
| ---------------------------- | ------------------------------------------------------------ |
| Ab 15.09.2025                | Ausgabe der Aufgabenstellung durch Dozierende an Studierende |
| Bis 19.12.2025 bis 18:00 Uhr | Abgabe der Projektarbeit in Complesis                        |
| 03.02.2026                   | Prov. Notenfreischaltung für Studierende ab 18:00 Uhr        |
| 27.02.2026                   | Def. Notenfreischaltung für Studierende ab 18:00 Uhr         |

### Arbeitsplatz

TE 617

## Hardware

- Hardware kann von Michael Peter und IMS bezogen werden
    - 2 x 3.5 Zoll Dronen
    - 2 x 5.0 Zoll Dronen
- Es können auch noch zwei weitere Dronen (1 x 3.5 Zoll, 1 x 5.0 Zoll) gebaut werden

## Ergebnisdarstellung

Ein Bericht im Umfang von 20-30 Seiten dokumentiert den Entwicklungsweg geeignet. Die Software wird in einem öffentlichen GitHub Repository abgelegt. Des weiteren soll eine Dokumentation mit Beispielen in Markdown erstellt werden. Die technische Dokumentation ist ebenso in Markdown zu verfassen.

## Fachgebiet

### Primäres Fachgebiet

Regelungstechnik und Advanced Control

### Weitere Fachgebiete

- Datenvisualisierung
- Digitale Signalverarbeitung
- Software Engineering

## Verschiedenes

### Studiengang

Systemtechnik

### Industriepartner

Keiner

### Institut / Zentren

Institut für Mechatronische Systeme (IMS)

### Interne Partner

Keine

## Informationslink

https://github.com/pichim/bf_controller_tuning

### Bemerkungen

- Grundsätzlich ist die Idee die PA aufbauend in einer BA zu erweitern

### Sprache

Deutsch / Englisch

## Studenten Detail

- Yuri Bianchi (biancyur)   - Automatiker   - Leader, Pilot, Regelungstechnik & Signalverarbeitung
- Dario Jurietti (juriedar) - Biolaborant   - Datenverarbeitung, Software
- Janick Dort (dortjan1)    - Elektroplaner - Datenverarbeitung, Pilot, Regelungstechnik & Signalverarbeitung

Die Studenten haben schon PM 1-4 zusammen absolviert und sind ein eingespieltes Team.
