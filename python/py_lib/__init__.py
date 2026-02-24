
"""
Betaflight Controller Tuning Analysis Library

This package provides tools for analyzing flight controller logs and tuning
PID controllers for quadcopters running Betaflight firmware.

Classes:
    FlightData: Load and preprocess flight log data
    GyroCtrlTuning: Controller analysis and tuning
    FlightAnalyzer: Spectral analysis (spectra, spectrograms)
    PlotUtils: Visualization tools

Modules:
    pidtuninglib: Helper functions for signal processing and control theory
"""

__version__ = '1.0.0'
__author__ = 'Janick Dort, Yuri Bianchi, Dario Jurietti'