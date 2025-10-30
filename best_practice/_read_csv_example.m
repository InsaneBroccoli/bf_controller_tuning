%% Best Practice Example: Reading CSV Files in MATLAB
% This example demonstrates how to read a CSV file and extract strings and numbers

%% Method 1: Using readtable() - RECOMMENDED for mixed data types
% This is the most straightforward approach for CSV files with both strings and numbers

% Read the CSV file
T = readtable('data.csv');

% View the table
disp(T);

% Extract specific columns
strings_column = T.StringColumn;  % Access by column name
numbers_column = T.NumericColumn;

% Convert to cell array if needed
string_data = table2cell(T(:, 'StringColumn'));
numeric_data = table2cell(T(:, 'NumericColumn'));

%% Method 2: Using readmatrix() with detectImportOptions() - For more control
% Better for handling specific data types and formats

opts = detectImportOptions('data.csv');

% Specify data types for each column
opts = setvartype(opts, 'StringColumn', 'string');
opts = setvartype(opts, 'NumericColumn', 'double');

% Read the file with options
T = readtable('data.csv', opts);

%% Method 3: Using csvread() - For numeric data only (legacy)
% NOTE: csvread is not recommended for mixed data; use readmatrix instead

data = readmatrix('data.csv');  % Modern alternative to csvread

%% Method 4: Manual parsing for complex scenarios
% Use when you need fine-grained control

fid = fopen('data.csv', 'r');
line_num = 0;
strings = {};
numbers = [];

while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        line_num = line_num + 1;
        % Skip header if present
        if line_num == 1
            continue;
        end
        
        % Split by comma
        parts = strsplit(line, ',');
        
        % Extract string (trim whitespace)
        strings{line_num - 1} = strtrim(parts{1});
        
        % Extract number
        numbers(line_num - 1) = str2double(strtrim(parts{2}));
    end
end
fclose(fid);

%% Practical Example with Sample Data
% Create a sample CSV file for demonstration

sample_data = {'Name,Age,Score', ...
               'Alice,25,95.5', ...
               'Bob,30,87.3', ...
               'Charlie,28,92.1'};

% Write sample data to CSV
fid = fopen('sample.csv', 'w');
for i = 1:length(sample_data)
    fprintf(fid, '%s\n', sample_data{i});
end
fclose(fid);

% Read the CSV
T = readtable('sample.csv');

% Extract data
names = T.Name;        % Strings
ages = T.Age;          % Numbers
scores = T.Score;      % Numbers

% Display results
fprintf('Names: ');
disp(names);
fprintf('Ages: ');
disp(ages);
fprintf('Scores: ');
disp(scores);

%% Best Practices Summary
% 1. Use readtable() for mixed data types - it handles strings and numbers automatically
% 2. Use readmatrix() for numeric-only data
% 3. Use detectImportOptions() when you need to specify data types
% 4. Always handle missing data with 'MissingRule' option
% 5. Use try-catch blocks for error handling in production code

%% Advanced: Handling Missing Data
opts = detectImportOptions('data.csv');
opts = setvaropts(opts, 'MissingRule', 'omitrow');  % Remove rows with missing data
T = readtable('data.csv', opts);

%% Example with Error Handling
try
    T = readtable('data.csv');
    disp('File read successfully');
catch ME
    fprintf('Error reading file: %s\n', ME.message);
end