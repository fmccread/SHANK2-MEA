% Setting up directories 
base_path = 'G:\Synchrony Analysis\CorSE Analysis 2023'; % folder with .raw files
addpath(genpath('G:\Synchrony Analysis\CorSE Analysis 2023\MATLAB functions')); % folder with custom matlab functions
cd(base_path)
mkdir("Analyzed Data")

%CorSE Parameters
fs = 12500;
WS = 1 * fs;
overlap = 0.5;

%weeks to analzye
weeks_analyze = [4,5,6,7,8];

%load plate information
plate_info = readtable('SHANK2_recording_info.csv', 'Delimiter', ',');

%get raw files
rawfiles = dir([base_path '\**\*.raw']);

%loop though each file
for fi = 1:length(rawfiles)
    
    %get file to load
    rawfile = [rawfiles(fi).folder '\' rawfiles(fi).name];

    %get plate and recording day
    splitName = split(rawfiles(fi).name, "_");
    plate = strrep(splitName{2},'-','');
    plate = str2num(strrep(plate,'Plate',''));
    day = str2num(strrep(splitName{3},'D',''));
    week = day2week(day);

    %skip recording if not in weeks being analyzed
    if ~ismember(week,weeks_analyze)
        disp("skipping file: Not in selected weeks to analyze")
        continue
    end

    %make new directory for plate data
    mkdir(strcat("Analyzed Data\Plate ", num2str(plate)))

    %get wells to analyze
    wells = string(table2array(plate_info(plate_info.Plate == plate & plate_info.Day == day & plate_info.Analyze == "YES", "Well")));

    %get line names
    lines = string(table2array(plate_info(plate_info.Plate == plate & plate_info.Day == day & plate_info.Analyze == "YES", "Line")));
    lines = strrep(lines,'19-2-2 WT','CTRL');
    lines = strrep(lines,'19-2-2 KO','KO');

    %indicies of wells to load
    [idx1, idx2] = well_indicies(wells);
    
    %load data
    disp(strcat('Analyzing recording', " ", num2str(fi), " of ", num2str(length(rawfiles)), ": ", rawfiles(fi).name))
    MEA_Data = AxisFile(rawfile).DataSets.LoadData(wells);

    % analyzing each well
    for ww = 1:length(wells)
        
        %metadata
        well = wells(ww);
        line = lines(ww);

        %get voltage vector
        VoltVec = [MEA_Data{idx1(ww),idx2(ww),:,:}];
        VoltVec = VoltVec.GetVoltageVector;

        %organizing into DataCell format 
        DataCell = cell(64,7);
        for k = 1:64
            DataCell{k,1} = k;
            DataCell{k,7} = VoltVec(:,k);
        end

        %display
        disp(strcat("Well ", num2str(ww)))
        
        %run CorSE analysis
        tic
        CorrData = par_CorSE_function_2023(DataCell, fs, WS, overlap);
        toc

        %saving data
        fname = strcat(base_path, "\Analyzed Data\Plate ", num2str(plate), '\Plate-', num2str(plate), '_D', num2str(day), "_", well, "_", line, '_CorrData.csv');
        csvwrite(fname,CorrData);

    end

end
