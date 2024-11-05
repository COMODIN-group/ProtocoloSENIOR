
function ProcessToEEGLabv3
    
    global ALLEEG;
    global EEG;
    global CURRENTSET;
    global ALLCOM;
    
    ALLEEG = [];
    EEG=[];
    CURRENTSET=[];
    ALLCOM=[];
        
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    output = {};
    nfiles = 0;
    listDir=dir(pwd);
    %Recorremos la lista de archivos
    listDir={listDir.name};
    for i=1:size(listDir,2)
            fileName=listDir(1,i);
            %solo procesamos los .dat y .gdf
            if ((~isempty(strfind(horzcat(fileName{:}), '.dat')))||(~isempty(strfind(horzcat(fileName{:}), '.gdf'))))
                output{end+1} = processFile(fileName);
                nfiles = nfiles + 1;
            end
    end
    labelrow = {'File'};
    labelrow{end+1} = {'% of signal kept'};
    labelrow{end+1} = {'ICs kept'};
    labelrow{end+1} = {'Interpolated channels'};
  
    output = reshape(output, [nfiles 1]);
    output = [{labelrow}; output];
    T = cell2table(output(1:end,:));
    % Write the table to a CSV file
    writetable(T,'EEG-Preprocessing-Statistics.csv','Delimiter','tab', 'WriteVariableNames',false);
    
    
end

function outputstr = processFile(filename)

        
    display(sprintf('File: %s ============================================================================>',horzcat(filename{:})));
        
    if (~isempty(strfind(horzcat(filename{:}), '.dat')))
        auxfilename = strrep(filename,'.dat','');
        EEG = pop_loadBCI2000(filename, 1);
    else
        auxfilename = strrep(filename,'.gdf','');
        filename = strcat('E:\Nacho\proyectos\IMDEA Alimentación\Buffer\',horzcat(filename{:}));
        %EEG = pop_biosig(filename, 'importevent','off', 'blockepoch','off');
        EEG = pop_biosig(filename);
    end
    datasetname = horzcat(auxfilename{:});
    
    outputstr = {datasetname};
    try
        % Select first 32 channels. The others are flat
        EEG = pop_select( EEG, 'nochannel',{'Channel 33','Channel 34','Channel 35','Channel 36','Channel 37','Channel 38','Channel 39','Channel 40','Channel 41','Channel 42','Channel 43','Channel 44','Channel 45','Channel 46','Channel 47','Channel 48','Channel 49','Channel 50','Channel 51','Channel 52','Channel 53','Channel 54','Channel 55','Channel 56','Channel 57','Channel 58','Channel 59','Channel 60','Channel 61','Channel 62','Channel 63','Channel 64'});
        EEG = eeg_checkset( EEG );
        % Channel locations
        EEG=pop_chanedit(EEG, 'lookup','C:\\Program Files\\MATLAB\\eeglab2022.0rc\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','load',{'E:\\Nacho\\proyectos\\IMDEA Alimentación\\Buffer\\channellocations.ced','filetype','autodetect'},'lookup','C:\\Program Files\\MATLAB\\eeglab2022.0rc\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
        EEG = eeg_checkset( EEG );
        % Resampling
        EEG = pop_resample(EEG, 256);
        EEG = eeg_checkset( EEG );
        % ASR
        beforeASRsize = EEG.pnts;
        EEG = eeg_checkset( EEG );
        EEG = clean_artifacts(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass',[0.25 0.75] ,'BurstCriterion',20,'WindowCriterion','off','BurstRejection','on','Distance','Euclidian');
        EEG = eeg_checkset( EEG );
        outputstr{end+1} = {(EEG.pnts/beforeASRsize) * 100};
        % Filter
        EEG = pop_eegfiltnew(EEG, 'locutoff',3,'hicutoff',31,'plotfreqz',0);
        EEG = eeg_checkset( EEG );
        % ICA
        EEG = pop_runica(EEG, 'icatype', 'runica');
        EEG = eeg_checkset( EEG );
        % ICLabel artifact removal
        EEG = pop_iclabel(EEG, 'default');
        EEG = eeg_checkset( EEG );
        EEG = pop_icflag(EEG, [NaN NaN;0.7 1;0.5 1;0.5 1;0.4 1;0.4 1;0.5 1]);
        EEG = eeg_checkset( EEG );
        EEG = pop_subcomp( EEG, [], 0, 0);
        outputstr{end+1} = size(EEG.etc.ic_classification.ICLabel.classifications, 1);
        % PREP
        EEG = pop_prepPipeline(EEG, struct('ignoreBoundaryEvents', true,'lineNoiseChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'lineFrequencies', [50  100], 'Fs', 256, 'p', 0.01, 'fScanBandWidth', 2, 'taperBandWidth', 2, 'taperWindowSize', 4, 'pad', 0, 'taperWindowStep', 1, 'fPassBand', [0  128], 'tau', 100, 'maximumIterations', 10, 'referenceChannels', [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32], 'evaluationChannels', [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32], 'rereferencedChannels', [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32], 'ransacOff', false, 'ransacSampleSize', 50, 'ransacChannelFraction', 0.25, 'ransacCorrelationThreshold', 0.75, 'ransacUnbrokenTime', 0.4, 'ransacWindowSeconds', 5, 'srate', 256, 'robustDeviationThreshold', 5, 'correlationWindowSeconds', 1, 'highFrequencyNoiseThreshold', 5, 'correlationThreshold', 0.2, 'badTimeThreshold', 0.01, 'maxReferenceIterations', 4, 'referenceType', 'Robust', 'reportingLevel', 'Verbose', 'interpolationOrder', 'Post-reference', 'meanEstimateType', 'Median', 'reportMode', 'skipReport', 'publishOn', false, 'sessionFilePath', '.\Report.pdf', 'summaryFilePath', '.\Summary.html', 'consoleFID', 1, 'cleanupReference', false, 'keepFiltered', true, 'removeInterpolatedChannels', false)); 
        % EEG = pop_prepPipeline(EEG, struct('lineNoiseChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'lineFrequencies', [50  100], 'Fs', 512, 'p', 0.01, 'fScanBandWidth', 2, 'taperBandWidth', 2, 'taperWindowSize', 4, 'pad', 0, 'taperWindowStep', 1, 'fPassBand', [0  256], 'tau', 100, 'maximumIterations', 10, 'reportMode', 'skipReport', 'publishOn', false, 'sessionFilePath', '.\Report.pdf', 'summaryFilePath', '.\Summary.html', 'consoleFID', 1, 'cleanupReference', false, 'keepFiltered', true, 'removeInterpolatedChannels', false, 'ignoreBoundaryEvents', false, 'detrendChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'detrendCutoff', 1, 'detrendStepSize', 0.02, 'detrendType', 'High Pass', 'referenceChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'evaluationChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'rereferencedChannels', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32], 'ransacOff', false, 'ransacSampleSize', 50, 'ransacChannelFraction', 0.25, 'ransacCorrelationThreshold', 0.75, 'ransacUnbrokenTime', 0.4, 'ransacWindowSeconds', 5, 'srate', 512, 'robustDeviationThreshold', 7, 'correlationWindowSeconds', 1, 'highFrequencyNoiseThreshold', 5, 'correlationThreshold', 0.3, 'badTimeThreshold', 0.01, 'maxReferenceIterations', 4, 'referenceType', 'Robust', 'reportingLevel', 'Verbose', 'interpolationOrder', 'Post-reference', 'meanEstimateType', 'Median', 'samples', 93120)); 
        EEG = eeg_checkset( EEG );
        % Average rereference
        EEG = pop_reref( EEG, []);
        EEG = eeg_checkset( EEG );
        %Channel interpolation by spectrum > 1 std
        EEG = eeg_checkset( EEG );
        [auxEEG, indelec, measure, com] = pop_rejchan( EEG, 'threshold', 1, 'measure', 'spec', 'norm', 'on', 'freqrange', [4 30]);
        outputstr{end+1} = size(union(EEG.etc.noiseDetection.interpolatedChannelNumbers, indelec), 2);
        EEG = eeg_checkset( EEG );
        EEG = pop_interp(EEG, indelec, 'spherical');
        EEG = eeg_checkset( EEG );
        % Processed file save
        EEG = pop_saveset( EEG, 'filename',strcat(datasetname,'_preprocesed.set'));
        % Processed file save for LORETA
        pop_export(EEG,strcat(datasetname,'_preprocesed-for-LORETA.txt'),'transpose','on','elec','off','time','off','precision',5);
        % Spectrum figure save
        figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', [4 7 10 15 21 28], 'freqrange',[4 30],'electrodes','off'); title(datasetname);
        saveas(gcf,strcat(datasetname,'.png'));
        close;
    catch
        %display(sprintf('ERROR: %s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',horzcat(filename{:})));
        outputstr = {datasetname};
        outputstr{end+1}={'Error!'}; 
        outputstr{end+1}={'Error!'}; 
        outputstr{end+1}={'Error!'}; 
    end
    
    EEG=[]; 
    eeglab redraw;
    outputstr = horzcat(outputstr{:});
    
end