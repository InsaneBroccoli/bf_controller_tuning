
# Betaflight Controller tuning
## Info
```
students:   Yuri Bianchi
            Janick Dort 
            Dario Jurietti

time period: september - december 2025
supervisor: Michael Peter, Ruprecht Altenburger
institution: ZHAW School of Engineering
```

---

### Was de Michi will
- fokus auf die arbeitsschritte damit verständlich ist was und wie wir es gemacht haben
    - wie wenn wir neu in einer firma anfangen mit einer aufgabe und ein vorarbeiter eine technische doku hinterlassen hat
- keine zu tiefen technischen details, wie herleitungen etc. 
    - dazu der Verweis auf GitHub repo für genauere Informationen
- Übung für BA -> unterteilte Abschnitte für individuelle Bewertung
- 

---

# Abstract
This project focuses on the further development of a method for offline tuning of rate controllers for FPV racing drones based on the open-source firmware Betaflight. The quality of rate control is crucial for flight stability, suppression of disturbances such as prop wash, and the agility of the drone. Until now, tuning has mostly been done through time-consuming trial and error.
The aim of the project is to develop a modular, open-source library for offline tuning based on an existing proof of concept. The library will be implemented in MATLAB and Python and supplemented by detailed Markdown documentation with examples on GitHub. The solution should offer the community a simple way to efficiently optimise drone parameters and thereby improve flight performance.
The results will be published in a public GitHub repository and documented in a report. In addition to software development, the project also includes practical work with FPV drones to validate the methods developed.

# Foreword

### Zweck und Motivation
- Warum wurde das Projekt durchgeführt?
- Welche Bedeutung hat es für dich oder das Team?
### Rahmenbedingungen
- In welchem Kontext entstand die Arbeit (Studium, Projektarbeit, Kooperation)?
- Eventuell kurze Erwähnung des Fachgebiets oder der Zielsetzung.
### Dank
- Danksagung an Betreuer, Unterstützer oder Personen, die geholfen haben.
### Persönliche Note
- Optional: Erfahrungen, Herausforderungen oder Lernaspekte.

# Table of contents



# 1. Introduction
As part of this project work…

## 1.1 Initial situation
The open-source firmware Betaflight is the most widely used platform for FPV racing drones and is continuously being developed by a large community. The quality of rate control is crucial for flight performance, as it ensures stability, agility and the suppression of disturbances such as prop wash.
Currently, most pilots tune these controller parameters primarily through trial and error, which is time-consuming and difficult for inexperienced pilots. A proof of concept (developed by Michael Peter) already exists that enables offline tuning of rate controllers based on a measurement data model (chirp signal experiment during flight). This concept is currently coded in MATLAB, which is not a widely used Software in the community and therefore rather inaccessible to the open world.

## 1.2 task and objective
The project aims to further develop the existing proof of concept, to make the tool more accessible to the Betaflight community and create an easy-to-use solution to improve the flight performance of FPV drones and replace the time-consuming trial-and-error process of tuning.
- has a clean, modularly expandable structure
- is provided with detailed Markdown how-to documentation and examples for users

# 2. Preparation
## 2.1 Chirp

### 2.1.1 Chirp concept
A chirp is an input signal that sweeps continuously over a defined range (e.g. from a low frequency to a high frequency). This signal excites the system on multiple frequencies in one run and thus lets you observe the frequency response and step response, the foundation of controller tuning.

Mathematically, a chirp signal can be represented as a sinusoidal wave with instantaneous frequency that varies linearly with time, such as $x(t) = \sin(2\pi (f_0 t + \frac{k}{2} t^2))$, where $(f_0)$
is the starting frequency and $(k)$ is the chirp rate.

### 2.1.2 Chirp on drone
In the context of tuning a drone's flight controller, an exponential sweep is used. This is due to the dominant low-frequency dynamics, where low frequencies exhibit more predictable responses, allowing the sweep to allocate more time to higher frequencies for better resolution of system behavior. Additionally, multirotor vibration resonances are typically spaced logarithmically, making exponential sweeps more suited for efficient excitation across the frequency range.

The process involves these steps:
  1. Signal Injection: The chirp signal is added to the pilot's stick input. This composite signal becomes the setpoint for the PID controller.
  2. System Response: The drone attempts to follow this rapidly changing setpoint. The motors react, and the drone will oscillate at varying frequencies corresponding to the chirp signal.
  3. Data Logging: While the chirp is active, high-resolution data from the gyroscope and the motor outputs are recorded using Betaflight's Blackbox feature. This data captures how the drone responds to the control inputs at different frequencies.

The recorded log file, containing the input chirp and the corresponding gyroscope output, is then used for offline analysis to determine the frequency response of the system.

### 2.1.3 Chirp programming
The implementation of the chirp signal in Betaflight is straightforward. You call the function chirpInit() and pass:
 - a struct of type chirp_t defined in chirp.h
 - starting frequency
 - end frequency
 - signal length
 - loop-time in microseconds

This initialises the values in the data struct and:
 - converts loop period from microseconds to seconds
 - calculates the number of samples in the signal duration
 - calculates the exponential growth factor
 - precomputes the constants used to get the phase from the frequency evolution
 - calls chirpReset() function

The chirpReset() function sets:
 - internal counter to 0
 - bool isFinished to false
then calls the function chirpResetSignals() to set signal-related fields to 0.

Once the chirp is initialised use the chirpUpdate() function once per loop iteration. This function checks:
 - if bool isFinished is true -> return false
 - else if internal counter is equal to the number of samples inside the chirps period -> set bool isFinished true, call chirpResetSignals(), return false
if none of the above is true the 

# 3. Theory

## 3.1 signal processing


## 3.2 control system


## 3.3 evaluation


# 4. Working steps and Problems
## 4.1 Simulations


## 4.2 Building an FPV drone


## 4.3 tuning an FPV drone
As part of understanding the workflow we tuned an FPV drone by ourselves. The goal was to adjust the PID parameters so that the step response and `disturbance rejection` could be analyzed using the proof-of-concept and later our restructured code.

Initially, we deliberately applied a poor tuning configuration to observe its effects. During the first field test with this setup, it quickly became apparent that the tuning was excessively poor. The drone began oscillating aggressively and climbed rapidly. As the drone was out of control, the kill switch had to be activated, resulting in a crash that broke one of the four arms. Fortunately, the damage was minor and repaired quickly.

After this incident, we decided not to repeat flights with intentionally bad tuning. Instead, we continued using the default Betaflight parameters as a baseline. `With these the drone flies reliably.` However, during demanding maneuvers such as flips or propwash manouvers, it becomes evident that there is room for improvement in terms of stability and responsiveness.

The tuning process was carried out on a larger 6-inch FPV drone provided by one of the team members. The drone was freshly flashed with Betaflight v12, meaning all filters and settings were at factory defaults. Achieving an optimal tune required several attempts due to challenges with filter adjustments and occasional logging issues. Nevertheless, the drone was successfully tuned, and the improvement was noticeable. During propwash manoeuvres, the drone remained significantly more stable, and flips were executed with greater precision at the stop.

Our tuning objective was to achieve a fast step response without overshoot, while maintaining good compliance and sensitivity. It is important to note that there is no universally “perfect” tune. Settings depend heavily on pilot preferences. Racing pilots for example maybe want more responsiveness for competitive performance, whereas freestyle pilots want smoothness and flight feel. 

Since none of our team members are professional pilots, we focused on achieving robustness and a step response without overshoot. This tuning approach is probably somewhere between what a racing pilot and a freestyle pilot would prefer. Robustness in this context is important for practical scenarios, such as when a propeller is slightly bent, the battery sits slightly off center or when an additional component like a camera is mounted on the drone.

## 4.4 Decode

## 4.5 Code Structure and programming style
From the start of the project, we recognized that a clear and maintainable code structure would be essential. The initial proof-of-concept code already contained a rough separation into topics such as data I/O and step response, but everything was packed into a single file, making it difficult to navigate. Our first step was to comment out the proof-of-concept code and reorganize it into thematic sections. This was the foundation for a more organized layoutand we came to an initial directory structure:

```
bf_controller_tuning/
│
├── main.m
├── data_io/
│   ├── load_data.m
│   ├── extract_header_information.m
│   └── preprocess_data.m
├── flight_parameters/
│   ├── get_pid_parameters.m
│   ├── scale_pid_parameters.m
│   └── print_parameters.m
├── spectrogram/
│   ├── calculate_spectrogram.m
│   ├── plot_spectrogram.m
│   └── show_spectrogram.m
├── spectra/
│   ├── estimate_spectra.m
│   ├── plot_spectra.m
│   └── show_spectra.m
├── frequency_response/
│   ├── estimate_frequency_response.m
│   ├── plot_bode.m
│   └── show_frequency_response.m
├── controller_analysis/
│   ├── calculate_closed_loop.m
│   ├── calculate_step_response.m
│   └── show_controller_analysis.m
└── utils/
    ├── apply_rotfiltfilt.m
    ├── downsample_frd.m
    └── get_my_colors.m
```
This draft was never intended to be final. Rather, it served as an initial attempt to bring clarity and separation between different functional areas.

### 4.5.1 OOP vs. Functional Approach
At this stage, we discussed whether to implement the entire project using object-oriented programming (OOP) or stick to a modular, function-based approach. Initially, we considered OOP as a challenge, since most of us were not very experienced with it. We started writing some components as classes, but quickly realized that only one team member really was comfortable with OOP, while the others struggled. This made us reconsider.
Our conclusion was to use OOP selectively—only for parts that benefit from reuse and encapsulation, such as data loading and plotting, while keeping the rest function-based for simplicity. This hybrid approach allowed us to maintain flexibility without overcomplicating the design.

After several attempts, we agreed on the following principles:

- A main script (main.m) for user input, figure selection, filter settings, and PID tuning
- A core class to handle all calculations (essentially the logic from the proof of concept, but without plotting)
- A plot utilities class to centralize all plotting functions for consistency

Existing library functions from the proof of concept were reused without modification. The thematic grouping included:

**Isch das eher d ordnerstruktur oder zellt das als Codestruktur??**
```
bf_controller_tuning/
│
├── class/
│     ├─ main_class.m
│     └── plot_utils.m
├── lib/
│     ├─ apply_rotfiltfilt.m
│     ├─ calculate_closed_loop.m
│     ├─ calculate_controllers.m
│     ├─ calculate_step_response_from_frd.m
│     ├─ calculate_transfer_functions.m
│     ├─ downsample_frd.m
│     ├─ estimate_frequency_response.m
│     ├─ estimate_spectra.m
│     ├─ estimate_spectrogram.m
│     ├─ expand_multiple_figure_nr.m
│     ├─ extract_header_information.m
│     ├─ get_chirp_signals.m
│     ├─ get_fcut_from_D_and_fcenter.m
│     ├─ get_fcut_from_exp.m
│     ├─ get_filter.m
│     ├─ get_ind_eval.m
│     ├─ get_my_colors.m
│     ├─ get_notch_Q.m
│     ├─ get_pid_scale.m
│     └─ get_switch_case_text_from_para.m
├── logs/
└── main.m
```
This structure worked well and kept the main script simple for the user. However, after a review meeting, our instructor suggested further modularization. The reason: we had created a “god class,” which limited flexibility. For example, if someone wanted to use the tool only for plotting spectrograms without any drone-specific logic, the current design made that difficult.

### 4.5.2 Final Structure Concept
To be written...


## 4.6 Python in MATLAB
The code has been divided into logically related sections in order to clearly separate topics such as spectral analysis or frequency response estimation and calculation. Dividing the code into thematic sections helps to improve readability and clearly present the individual topics. As we thaught the MATLAB code was structured good enough, we wanted to start with the conversion to python. First we considered starting from scratch and rewrite the whole code, but it was made clear, that this is hard to test within the rewriting process. The idea then was to convert all the MATLAB functions in the `\lib` to python first and then call these rewritten python functions from the `gyro_ctrl_tuning` class in MATLAB. The advantage of doing so, gave us the ablity to test every function, weather they were converted correctly or not. In a second step, the `gyro_ctrl_tuning` and lastly the `plot_utils` class would have been converted to python. Unfortunately we ran into the problem, $that in python no FRD objects exist$. This led us to write helper functions, which converted the MATLAB FRD objects into NumPy-arrays using the numpy library. One array holding the frequencies and one the response data. The same process had to be done in the opposite direction, from python to MATLAB. It quickly became clear, that this was too much effort. This is because later in python, these helperfunctions are not necessary anymore, as we would hold this data in two arrays anyway. 

As Dario started an easy python version of the tool, the plan was to just build up on this. 

# 5. Results

# 6. Discussion and next steps
## 6.1 Discussion

## 6.2 Next steps
- Further develop the tool as part of a bachelor’s thesis
- Convert the code to python

## 6.3 Conclusion


# 7. Directories
