% csv_surf_converter_script.m
% Select input CSV and output SURF file via dialog boxes, then convert

% Prompt user to select the input CSV file
[input_file, input_path] = uigetfile('*.csv', 'Select the input CSV file');
if isequal(input_file,0)
    disp('User canceled input file selection.');
    return
end
input_csv = fullfile(input_path, input_file);

% Prompt user to specify the output SURF file name
[output_file, output_path] = uiputfile('*.surf', 'Save output SURF file as');
if isequal(output_file,0)
    disp('User canceled output file selection.');
    return
end
output_surf = fullfile(output_path, output_file);

% --- Now process the CSV and write SURF ---
try
    % Read the whole file as text lines
    fid = fopen(input_csv, 'r');
    if fid == -1
        error('Cannot open input file.');
    end
    lines = {};
    tline = fgetl(fid);
    while ischar(tline)
        lines{end+1} = tline; %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    % Find the line index for "Airfoil surface,"
    idx = find(strcmp(strtrim(lines), 'Airfoil surface,'), 1);
    if isempty(idx)
        error('Airfoil surface section not found.');
    end

    % Data start two lines after 'Airfoil surface,' (skip header line)
    data_start = idx + 2;

    % Read until blank line or next section (line without a comma)
    coords = [];
    for i = data_start:length(lines)
        line = strtrim(lines{i});
        if isempty(line) || ~contains(line, ',')
            break;
        end
        parts = strsplit(line, ',');
        x = str2double(parts{1});
        y = str2double(parts{2});
        coords(end+1, :) = [x, y]; %#ok<AGROW>
    end

    % Normalize X and Y by chord length (100 mm)
    chord = 100;
    coords(:,1) = coords(:,1) / chord;
    coords(:,2) = coords(:,2) / chord;

    % Write to .surf file with formatting
    fid_out = fopen(output_surf, 'w');
    if fid_out == -1
        error('Cannot open output file.');
    end
    for i = 1:size(coords, 1)
        fprintf(fid_out, '%10.6f %10.6f \n', coords(i,1), coords(i,2));
    end
    fclose(fid_out);

    fprintf('Converted %s to %s\n', input_csv, output_surf);
catch ME
    disp('Error during conversion:');
    disp(ME.message);
end
