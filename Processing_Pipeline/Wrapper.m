%% PREPROCESSING AND ANALYSIS WRAPPER
%clear and close everything
ccc

%%%%
% description of the dataset settings
% 1) NewBaseline: this has the window of 1024ms, and the baseline is 712-512ms. Epoched to the entrainers
% 2) NewBaseline2: this has the window of 512ms, and the baseline is 456-256ms. Epoched to the entrainers
% 3) Target: Windowsize of 512ms, no baseline epoched to target onset and divided into Valid and Invalid targets
% 4) TargetERP: no time-f, locked to first erp, filter at 50 Hz
%%%%


%% Settings for loading the raw data
%Datafiles must be in the format exp_participant, e.g. EEGexp_001.vhdr
exp.name = 'MultiFreq';
% List of participants' ids
exp.participants = {'001';'002';'003';'004';'005';'008';'009';'010';'011';'012';'014';'015';'016';'017';'018';'019';'023';'024';'025';'026';'027';'028';'029';'030';'031';'032';'033';'034';'035';'036';'037';'038';'039';'040';'041';'042'};
exp.pathname = 'M:\Data\MultiFreq\EEG'; %path of EEG data
% The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
exp.settings = 'Baseline';

%% Choose how to organize your datasets. 
% The sets are different trial types that contain comparable events. 
% Examples of different sets: 10Hz/12Hz stimulation on different trials; emotional vs non-emotional images on different trials.
% Each row is a different set. 
exp.setname = {'Baseline'}; %name the rows

exp.conds = ''; %conds is usually used for comparing the same type of trials
                %under different conditions (e.g., stimulation vs sham)

%% Events and event labels
%Events are what happen within each trial
%Examples of different events. Missing vs present targets, valid vs invalid targets,
%detected vs undetected responses.

%%Each column is a different event. You can collect multiple triggers into one
%%event with square brackets [].

%%%For each condition (lag 1-4 in this case), numbers correspond to
%%%triggers that will be kept for each condition. All other triggers will
%%%be removed
exp.events = {[11 21], [13 23]};    %must be matrix (sets x events)       
exp.event_names = {'Eyes Closed','Eyes Open'}; %name the columns

% Each item in the exp.events matrix will become a seperate dataset, including only those epochs referenced by the events in that item. 
%e.g. 3 rows x 4 columns == 12 datasets/participant

%% Blink Correction
% the Blink Correction wants dissimilar events (different erps) seperated by 
% commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
exp.selection_cards = {'11 21','13 23'}; %must be list == length(exp.setname)

%% Artifact rejection. 
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold 
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
% exp.postocularthresh = [ ]; %Second happens after. Leave blank [] to skip
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip

%% Electrode location
%Where are your electrodes? (.ced file)
exp.electrode_locs = 'M:\Analysis\OrientWheel\EOG-electrode-locs-32_orientwheel.ced';
% electrode information
exp.electrode = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
exp.elec_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

%% Re-referencing the data
exp.refelec = 1; %which electrode do you want to re-reference to?
exp.brainelecs = [2:32]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Filter the data?
exp.filter = 'on'; %filter all files
exp.hicutoff = 50; %higher edge of the frequency pass band (Hz)
exp.locutoff = 0.1; %lower edge of the frequency pass band (Hz)

%% FFT/Wavelet settings
% How long is your window going to be? (Longer window == BETTER frequency 
% resolution & WORSE time resolution)
exp.winsize = 512; %use numbers that are 2^x, e.g. 2^10 == 1024ms

% Baseline will be subtracted from the power variable. It is relative to 
% your window size. Can use just NaN for no baseline
%e.g., [-200 0] will use [-200-exp.winsize/2 0-exp.winsize/2]; 
exp.erspbaseline = NaN;
% exp.erspbaseline = [-400 -200];

% Instead of choosing a windowsize, you can choose a number of cycles per 
% frequency for standard wavelet analysis: usually [3] for better temporal
% precision or [6] for better frequency precision.
% If [wavecycles factor], wavelet cycles increase with frequency beginning 
% at wavecyles. See "help popnewtimef"
% exp.cycles = [0]; %leave it at 0 to use FFT
exp.cycles = [2 0.8]; %number of cycles 

% Choose number of output times
exp.timesout = 300; %200 is usually used

% Set sampling factor for frequencies. 
% when exp.cycles==0, frequency spacing is (low_freq/padratio). For wavelet,
% multiplies the # of output freqs by dividing their spacing (2 is default).
% higher values give smooth looking plot but at computational cost (16 is
% very high)
exp.padratio = 4;

% What frequencies to consider?
% exp.freqrange = [1 40]; 
exp.freqrange = [exp.cycles(1) 40]; %when doing wavelet

%% Epoching the data
exp.epoch = 'on'; %on to epoch data; off to load previous data
%%%indicates where you want to center your data (where time zero is)
exp.epochs = {}; %must be list == length(exp.setname)
exp.epochs_name = {};
exp.epochslims = [-1.5 1.5]; %in seconds; epoched trigger is 0 e.g. [-1 2]
exp.epochbaseline = [-200 0]; %remove the baseline for each epoched set, in ms. e.g. [-200 0] 


%% Time-Frequency settings
%Do you want to run time-frequency analyses? (on/off)
exp.tf = 'on';
%Do you want to use all the electrodes or just a few? Leave blank [] for 
% all (will use same as exp.brainelecs)
exp.tfelecs = [];

%Do you want to save the single-trial data? (on/off) (Memory intensive!!!)
exp.singletrials = 'on';
%Saving the single trial data is memory intensive. Just use the electrodes
% you need. 
exp.singletrialselecs = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 20 21 22 25 26];
exp.singtrlelec_name = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';...
    'P7';'P8';'P5';'P6';'P3';'P4';'CP1';'CP2';'C3';'C4';'FC1';'FC2'};

%% Save your pipeline settings
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. This will help you remember what settings were used in each dataset



%% Run Preprocessing
Preprocessing(exp) %comment out if you're only doing analysis
% Preprocessing_Example2(exp)

%% Run Analysis
%Don't want to change all the above settings? Load the settings from the saved .mat file.

%choose the data types to load into memory (on/off)
anal.segments = 'on'; %load the EEG segments?
anal.tf = 'on'; %load the time-frequency data?

anal.singletrials = 'on'; %load the single trial data?
anal.entrainer_freqs = [20; 15; 12; 8.5; 4]; %Single trial data is loaded at the event time, and at the chosen frequency. 

anal.tfelecs = []; %load all the electodes, or just a few? Leave blank [] for all.
anal.singletrialselecs = [2 3 4 6];

Analysis(exp,anal) % The Analysis primarily loads the processed data. It will attempt to make some figures, but most analysis will need to be done in seperate scripts.































