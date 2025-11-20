
# Betaflight Controller tuning
students: Yuri Bianchi, Janick Dort, Dario Jurietti

date: HS 2025

institution: ZHAW School of Engineering

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

Mathematically, a chirp signal can be represented as a sinusoidal wave with instantaneous frequency that varies linearly with time, such as 
\( x(t) = \sin(2\pi (f_0 t + \frac{k}{2} t^2)) \), where \( f_0 \)
is the starting frequency and $(k)$ is the chirp rate. In controller tuning, chirps are particularly useful for identifying system dynamics in real-time, such as in flight controllers. However, they can introduce nonlinear effects if the sweep is too rapid, potentially distorting the response."
### 2.1.2 Chirp on drone
### 2.1.3 Chirp programming



# 3. Theory


# 4. Working steps and Problems
## 4.1 Simulations

## 4.2 Code Structure

## 4.3 Building FPV drone

## 4.4 Trying to implement the conversion from log to csv
Get same Headerfile as with the bf blackbox viewer

# 5. Results

# 6. Discussion and next steps
## 6.1 Discussion

## 6.2 Next steps
Further develop the tool as part of a bachelor’s thesis
Convert the code to python

## 6.3 Conclusion


# 7. Directories
