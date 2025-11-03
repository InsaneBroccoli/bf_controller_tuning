function [strings, numbers] = extract_header(file_path)
  
  fid = fopen(file_path);
  line_num = 0;
  strings = {};
  numbers = [];

  while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
      line_num = line_num + 1;
      % Skip to line 11
      if line_num < 12
        continue;
      end

      % split by comma
      parts = strsplit(line, ',');

      % extract string (trim whitespace)
      strings{line_num - 1} = strtrim(parts{1});

      % extract number
      numbers(line_num - 1) = str2double(strtrim(parts{2}));
    end
  end
  fclose(fid);
