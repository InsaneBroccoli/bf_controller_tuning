clc; clear all; close all;
addpath ./my_lib/

folder_name = "testing";
file_name = "testing.csv";

file_path = fullfile(folder_name, file_name)
[names, numbers] = extract_header(file_path);