classdef extract_header_information
    properties (SetAccess = private)
        filePath (1,1) string
        para struct = struct()
        ind  struct = struct()
        Nheader (1,1) double = 0
    end

    methods
        % Check if a filepath is given
        function obj = extract_header_information(filePath)
            % nargin stands for number of inputs and is provided by MATLAB
            if nargin > 0
                obj.filePath = string(filePath);

                if ~isfile(obj.filePath)    % In case the file is not found
                    error("The file '%s' was not found.", obj.filePath);     % Error message
                end
            end
        end

        function obj = parseHeader(obj)
            % Open file for reading and treat file as text
            fid = fopen(obj.filePath, 'rt');
            if fid == -1
                error('The file could not be opened.');
            end
            cleaner = onCleanup(@() fclose(fid)); % Registers cleanup when the function finishes

            % Reset data for multiple runs
            obj.Nheader = 0;
            obj.para = struct();
            obj.ind  = struct();

            tline = fgetl(fid);                 % Get data from the first line of the file
            do_read_para   = false;             %
            lastHeaderLine = '';                % When we get no data from a line, the program sees this as the end of the file

            while ischar(tline)                 % As long as we get a char this loop will go on
                obj.Nheader = obj.Nheader + 1;  % Counter for the number of lines we get

                % Ignore every line until the line frameInterval and then start getting parameter data
                if ~isempty(regexp(tline, 'frameIntervalI', 'once'))
                    do_read_para = true;
                end

                % Stop here because after this you get the flight data
                if ~isempty(regexp(tline, 'loopIteration', 'once'))
                    lastHeaderLine = tline;      % Get the names of the different flight data
                    break;
                end

                if do_read_para
                    idx = regexp(tline, '",');   % Split the parameter name and data
                    if ~isempty(idx)
                        para_name  = tline(2:idx-1);   % Get the name of the parameter
                        para_value = tline(idx+2:end); % Get the parameter value
                        if ~isempty(para_value) && para_value(1) == '"'   % If the value is quoted
                            % Set the parameter name and the data together
                            try % 'magPID' '"40,,"' because this parameter is a special case
                                eval(['obj.para.(para_name) = [', para_value(2:end-1), '];']);
                            end
                        else
                            eval(['obj.para.(para_name) = [', para_value, '];']);
                        end
                    end
                end

                tline = fgetl(fid); % Go to the next line
            end

            % Parse the last scanned line (the "loopIteration,..." header with column names)
            if ~ischar(tline)
                error('Header parsing failed: reached end of file before finding the loopIteration line.');
            end

            idx = regexp(tline, ',');   % Get positions of commas in the last scanned line
            if isempty(idx)
                error('Header parsing failed: no commas found in the loopIteration line.');
            end

            ind_cntr = 1;       % Index counter

            % Handle first element, should be "loopIteration"
            ind_name = tline(2:idx(1)-2); %#ok
            eval(['obj.ind.(ind_name) = [', num2str(ind_cntr), '];']);

            % Handle all elements between
            for i = 1:length(idx)
                ind_cntr = ind_cntr + 1;
                if i < length(idx)
                    ind_name = tline(idx(i)+2:idx(i+1)-2);  % Get the name of every data field in loopIteration
                else
                    ind_name = tline(idx(end)+2:end-1);     % For the last name
                end

                % Get names of flight data
                if strncmp(ind_name, 'eRPM', 4)                 % keep eRPM special case
                    eval(['obj.ind.(ind_name(1:4))(str2double((ind_name(end-1))) + 1) = [', num2str(ind_cntr), '];']);
                elseif endsWith(ind_name, ']')
                    eval(['obj.ind.(ind_name(1:end-3))(str2double((ind_name(end-1))) + 1) = [', num2str(ind_cntr), '];']);
                else
                    eval(['obj.ind.(ind_name) = [', num2str(ind_cntr), '];']);
                end
            end
        end
    end
end
