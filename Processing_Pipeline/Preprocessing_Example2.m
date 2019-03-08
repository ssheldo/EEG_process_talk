function Preprocessing_Example2(exp)

% Not set-up to process data from the staircasing task

try
%     if parpool('poolsize') == 0
%         parpool OPEN 3;
%         parpool(3)
%     end
    
    nparts = length(exp.participants);
    nsets = length(exp.setname);
    
    % Replicating event names when exp.events is a matrix
%     if any(size(exp.event_names) ~= size(exp.events))
%         repfactor = int8(size(exp.events)./size(exp.event_names));
%         exp.event_names = repmat(exp.event_names, repfactor);
%     end

    % Replicating event triggers when exp.events is a matrix
    if isempty(exp.epochs) == 1
    %     exp.epochs = exp.events;
        exp.epochs = cellstr(num2str(cell2mat(reshape(exp.events,1,size(exp.events,1)*size(exp.events,2)) )'))';
    else
        exp.events = exp.epochs;
    end

    % Is epoch names not specified, use event names
    if isempty(exp.epochs_name) == 1
        exp.epochs_name = exp.event_names;
    else
        exp.event_names = exp.epochs_name;
    end
    
    
    for i_set = 1:nsets
        
        sprintf(exp.setname{i_set})
        
        % if folder doesn't exist yet, create one
        if ~exist([exp.pathname  '\' exp.setname{i_set} '\Segments\'])
            mkdir([exp.pathname  '\' exp.setname{i_set} '\Segments\']);
        end
        
        %initialize EEGLAB
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %subject numbers to analyze
        
        nevents = length(exp.events(i_set,:));
        
        %% Load data and channel locations
        for i_part = 1:nparts
%         for i_part = 12    
            
            part_name = ['subj' exp.participants{i_part}]; %this is cuz subject ids are in the form 'subj1' rather than '001' in orientation wheel
            
            if strcmp('on',exp.epoch) == 1
                
                sprintf(['Participant ' num2str(exp.participants{i_part})])
                
                %% Load a data file
                if strcmpi(part_name,'subj13') %subject's EEG file has different name
                    EEG = pop_loadbv(exp.pathname, [part_name '_staircase.vhdr']); %orientation wheel data only
                else
                    EEG = pop_loadbv(exp.pathname, [part_name '_orient.vhdr']); %orientation wheel data only
                end
                
                %% Load channel information
                EEG = pop_chanedit(EEG, 'load',{exp.electrode_locs 'filetype' 'autodetect'});
               
                %% Arithmetically re-reference to linked mastoid (M1 + M2)/2
                % only for brain electrodes (exclude mastoids & EOGs)
                for ii = exp.brainelecs(1):length(exp.brainelecs)
                    EEG.data(ii,:) = (EEG.data(ii,:)-((EEG.data(exp.refelec,:))*.5));
                end
                clear ii
                
                %% Filter the data
                if strcmpi(exp.filter,'on')
%                    EEG = pop_eegfilt( EEG, 0, 30, [], 0); %with low pass of 30
                   EEG = pop_eegfiltnew(EEG, exp.locutoff, exp.hicutoff); % filter function
                end
                
                %% Change markers so they can be used by the gratton_emcp script
                allevents = length(EEG.event);
                for i_event = 2:allevents %skip the first
                    EEG.event(i_event).type = num2str(str2num(EEG.event(i_event).type(2:end)));
                end

                %% The triggers are early
                [EEG] = VpixxEarlyTriggerFix(EEG);
                
                %% Extract epochs of data time locked to event
                %Extract data time locked to targets and remove all other events
                EEG = pop_epoch(EEG, exp.epochs, exp.epochslims, 'newname', [part_name '_epochs'], 'epochinfo', 'yes');
                %subtract baseline
                EEG = pop_rmbase(EEG, exp.epochbaseline);
  
                %% Get behavior data and add to EEG structure
                [error_deg] = getBEHdata_OrientWheel(part_name);
                % add error deg to epoch structure
                if length(EEG.epoch) == length(error_deg)... %make sure right BEH file
                        || strcmpi(part_name,'subj6') || strcmpi(part_name,'subj8') %or subject 6 and 8
                   EEG.error_deg = error_deg;
                end
                
                %% Reject practice trials from data
                % only specific to subject 8 % only when using these specific settings
                if (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byTargets_v3'))...
                    || (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byTargets_v1'))...
                    || (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byTargets_v2'))...
                    || (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'byTargets_v1'))...
                    || (strcmpi(part_name,'subj8') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))
                    rej_practice = zeros(1,length(EEG.epoch));
                    rej_practice(1,1:18) = 1; %mark the first 18 trials for removal
                    EEG = pop_rejepoch(EEG, rej_practice, 0);
                end
                clear rej_practice
                

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %% Artifact Rejection, EMCP Correction, then 2nd Rejection
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 1, trials with range >exp.preocularthresh uV
                if isempty(exp.preocularthresh) == 0
                    rejtrial = struct([]);
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],exp.preocularthresh(1),exp.preocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    % `````````````````````````````````````````````````````  
                    % additional rejection of trials for subject 7 due to technical issues (reviewed visually) 
                    if strcmpi(part_name,'subj7') && strcmpi(exp.settings,'filt_byTargets_v1') %only when using this specific settings 
                        EEG.reject.rejthresh(50)=1;
                        EEG.reject.rejthresh(65:71)=1;
                        EEG.reject.rejthresh(74:76)=1;
%                     elseif strcmpi(part_name,'subj7') && strcmpi(exp.settings,'filt_byTargets_v2') %only when using this specific settings 
%                         EEG.reject.rejthresh(50)=1;
%                         EEG.reject.rejthresh(65:67)=1;
%                         EEG.reject.rejthresh(70:72)=1;
                    end
                    % ````````````````````````````````````````````````````` 
                    rejtrial(i_set, 1).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % EMCP occular correction
                temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
                EEG = gratton_emcp(EEG, exp.selection_cards, {'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
                EEG.emcp.table %this prints out the regression coefficients
                EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Baseline again since EMCP changed it
                EEG = pop_rmbase(EEG,exp.epochbaseline);
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 2, trials with range >exp.postocularthresh uV
                if isempty(exp.postocularthresh) == 0
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],exp.postocularthresh(1),exp.postocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    rejtrial(i_set,2).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
                %% Additional rejection of trials for subject 7 due to technical issues (reviewed visually) 
                if (strcmpi(part_name,'subj3') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(16)=1;
                    EEG.reject.rejthresh(100)=1;
                    EEG.reject.rejthresh(167)=1;
                    EEG.reject.rejthresh(169)=1;
                    EEG.reject.rejthresh(173)=1;
                    EEG.reject.rejthresh(234)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj5') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(52)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj5') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(172)=1;
                    EEG.reject.rejthresh(185)=1;
                    EEG.reject.rejthresh(194)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
                elseif (strcmpi(part_name,'subj6') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(19)=1;
                    EEG.reject.rejthresh(58)=1;
                    EEG.reject.rejthresh(61)=1;
                    EEG.reject.rejthresh(68)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj6') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(70)=1;
                    EEG.reject.rejthresh(83)=1;
                    EEG.reject.rejthresh(113:115)=1;
                    EEG.reject.rejthresh(125)=1;
                    EEG.reject.rejthresh(136)=1;
                    EEG.reject.rejthresh(163)=1;
                    EEG.reject.rejthresh(191)=1;
                    EEG.reject.rejthresh(195)=1;
                    EEG.reject.rejthresh(201)=1;
                    EEG.reject.rejthresh(207:208)=1;
                    EEG.reject.rejthresh(218)=1;
                    EEG.reject.rejthresh(220:221)=1;
                    EEG.reject.rejthresh(238)=1;
                    EEG.reject.rejthresh(241)=1;
                    EEG.reject.rejthresh(244)=1;
                    EEG.reject.rejthresh(261:263)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
                elseif (strcmpi(part_name,'subj7') && strcmpi(exp.settings,'filt_byTargets_v2'))...
                    || (strcmpi(part_name,'subj7') && strcmpi(exp.settings,'filt_byTargets_v3')) %only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(39)=1;
                    EEG.reject.rejthresh(50:51)=1;
                    EEG.reject.rejthresh(64:72)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj9') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(152)=1;
                    EEG.reject.rejthresh(155)=1;
                    EEG.reject.rejthresh(234)=1;
                    EEG.reject.rejthresh(238)=1;
                    EEG.reject.rejthresh(261)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
                elseif (strcmpi(part_name,'subj10') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(48)=1;
                    EEG.reject.rejthresh(67)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                 elseif (strcmpi(part_name,'subj10') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(30)=1;
                    EEG.reject.rejthresh(79)=1;
                    EEG.reject.rejthresh(191:193)=1;
                    EEG.reject.rejthresh(250)=1;
                    EEG.reject.rejthresh(265)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
                 elseif (strcmpi(part_name,'subj11') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(48:49)=1;
                    EEG.reject.rejthresh(115)=1;
                    EEG.reject.rejthresh(121)=1;
                    EEG.reject.rejthresh(182)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj12') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(47)=1;
                    EEG.reject.rejthresh(68)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj12') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(48:49)=1;
                    EEG.reject.rejthresh(64)=1;
                    EEG.reject.rejthresh(96)=1;
                    EEG.reject.rejthresh(108)=1;
                    EEG.reject.rejthresh(118)=1;
                    EEG.reject.rejthresh(246)=1;
                    EEG.reject.rejthresh(250:252)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj13') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(33)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj13') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(10:11)=1;
                    EEG.reject.rejthresh(30)=1;
                    EEG.reject.rejthresh(151:152)=1;
                    EEG.reject.rejthresh(192)=1;
                    EEG.reject.rejthresh(230:231)=1;
                    EEG.reject.rejthresh(234)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
                elseif (strcmpi(part_name,'subj14') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(51)=1;
                    EEG.reject.rejthresh(55:56)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj14') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(91)=1;
                    EEG.reject.rejthresh(174)=1;
                    EEG.reject.rejthresh(189:190)=1;
                    EEG.reject.rejthresh(201)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj15') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(57)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                 elseif (strcmpi(part_name,'subj15') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(58:59)=1;
                    EEG.reject.rejthresh(81:83)=1;
                    EEG.reject.rejthresh(104)=1;
                    EEG.reject.rejthresh(149)=1;
                    EEG.reject.rejthresh(151)=1;
                    EEG.reject.rejthresh(153:154)=1;
                    EEG.reject.rejthresh(160:161)=1;
                    EEG.reject.rejthresh(197)=1;
                    EEG.reject.rejthresh(200:203)=1;
                    EEG.reject.rejthresh(205:206)=1;
                    EEG.reject.rejthresh(228)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj16') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(43)=1;
                    EEG.reject.rejthresh(53)=1;
                    EEG.reject.rejthresh(120)=1;
                    EEG.reject.rejthresh(230)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj17') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(55)=1;
                    EEG.reject.rejthresh(63)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj17') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(2)=1;
                    EEG.reject.rejthresh(47)=1;
                    EEG.reject.rejthresh(66)=1;
                    EEG.reject.rejthresh(101)=1;
                    EEG.reject.rejthresh(152)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj18') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(70)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj18') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(101)=1;
                    EEG.reject.rejthresh(257)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byTargets_v2'))...
                        || (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(31)=1;
                    EEG.reject.rejthresh(82:85)=1;
                    EEG.reject.rejthresh(93)=1;
                    EEG.reject.rejthresh(104:108)=1;
                    EEG.reject.rejthresh(115:116)=1;
                    EEG.reject.rejthresh(139)=1;
                    EEG.reject.rejthresh(146)=1;
                    EEG.reject.rejthresh(154:158)=1;
                    EEG.reject.rejthresh(181)=1;
                    EEG.reject.rejthresh(186:188)=1;
                    EEG.reject.rejthresh(191)=1;
                    EEG.reject.rejthresh(195:199)=1;
                    EEG.reject.rejthresh(220:221)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))...
                       || (strcmpi(part_name,'subj19') && strcmpi(exp.settings,'filt_byCatchTargets_v1')) %only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(14)=1;
                    EEG.reject.rejthresh(20)=1;
                    EEG.reject.rejthresh(25:26)=1;
                    EEG.reject.rejthresh(28)=1;
                    EEG.reject.rejthresh(41:42)=1;
                    EEG.reject.rejthresh(46)=1;
                    EEG.reject.rejthresh(58)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
                elseif (strcmpi(part_name,'subj21') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(47)=1;
                    EEG.reject.rejthresh(69)=1;
                    EEG.reject.rejthresh(71)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj21') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(30:31)=1;
                    EEG.reject.rejthresh(116)=1;
                    EEG.reject.rejthresh(118)=1;
                    EEG.reject.rejthresh(145:147)=1;
                    EEG.reject.rejthresh(175)=1;
                    EEG.reject.rejthresh(179:181)=1;
                    EEG.reject.rejthresh(192)=1;
                    EEG.reject.rejthresh(206:207)=1;
                    EEG.reject.rejthresh(212:213)=1;
                    EEG.reject.rejthresh(231)=1;
                    EEG.reject.rejthresh(239)=1;
                    EEG.reject.rejthresh(249)=1;
                    EEG.reject.rejthresh(251)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj22') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(4)=1;
                    EEG.reject.rejthresh(40)=1;
                    EEG.reject.rejthresh(95)=1;
                    EEG.reject.rejthresh(99)=1;
                    EEG.reject.rejthresh(104)=1;
                    EEG.reject.rejthresh(153)=1;
                    EEG.reject.rejthresh(245)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))...
                        || (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byCatchTargets_v1')) %only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(9)=1;
                    EEG.reject.rejthresh(18)=1;
                    EEG.reject.rejthresh(24)=1;
                    EEG.reject.rejthresh(33)=1;
                    EEG.reject.rejthresh(40)=1;
                    EEG.reject.rejthresh(60)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
                elseif (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byTargets_v2'))...
                       || (strcmpi(part_name,'subj23') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    %subject moved throughout experiment
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(10)=1;
                    EEG.reject.rejthresh(23)=1;
                    EEG.reject.rejthresh(26)=1;
                    EEG.reject.rejthresh(39)=1;
                    EEG.reject.rejthresh(42)=1;
                    EEG.reject.rejthresh(51)=1;
                    EEG.reject.rejthresh(57)=1;
                    EEG.reject.rejthresh(72:73)=1;
                    EEG.reject.rejthresh(87)=1;
                    EEG.reject.rejthresh(107)=1;
                    EEG.reject.rejthresh(110)=1;
                    EEG.reject.rejthresh(117)=1;
                    EEG.reject.rejthresh(132:133)=1;
                    EEG.reject.rejthresh(141)=1;
                    EEG.reject.rejthresh(149)=1;
                    EEG.reject.rejthresh(156)=1;
                    EEG.reject.rejthresh(158)=1;
                    EEG.reject.rejthresh(167)=1;
                    EEG.reject.rejthresh(170:172)=1;
                    EEG.reject.rejthresh(174)=1;
                    EEG.reject.rejthresh(178)=1;
                    EEG.reject.rejthresh(180:182)=1;
                    EEG.reject.rejthresh(188:189)=1;
                    EEG.reject.rejthresh(194)=1;
                    EEG.reject.rejthresh(201)=1;
                    EEG.reject.rejthresh(214)=1;
                    EEG.reject.rejthresh(220)=1;
                    EEG.reject.rejthresh(228)=1;
                    EEG.reject.rejthresh(231)=1;
                    EEG.reject.rejthresh(241)=1;
                    EEG.reject.rejthresh(246)=1;
                    EEG.reject.rejthresh(252:253)=1;
                    EEG.reject.rejthresh(258)=1;
                    EEG.reject.rejthresh(263:264)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
                elseif (strcmpi(part_name,'subj24') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(69)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);  
                 elseif (strcmpi(part_name,'subj24') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(119)=1;
                    EEG.reject.rejthresh(192)=1;
                    EEG.reject.rejthresh(228)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                 elseif (strcmpi(part_name,'subj25') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    %subject wore a lot of makeup
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(13)=1;
                    EEG.reject.rejthresh(17)=1;
                    EEG.reject.rejthresh(19)=1;
                    EEG.reject.rejthresh(70)=1;
                    EEG.reject.rejthresh(98)=1;
                    EEG.reject.rejthresh(101)=1;
                    EEG.reject.rejthresh(208)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);    
                elseif (strcmpi(part_name,'subj26') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(14:15)=1;
                    EEG.reject.rejthresh(53)=1;
                    EEG.reject.rejthresh(115)=1;
                    EEG.reject.rejthresh(240)=1;
                    EEG.reject.rejthresh(269)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);   
                elseif (strcmpi(part_name,'subj27') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(2)=1;
                    EEG.reject.rejthresh(7)=1;
                    EEG.reject.rejthresh(75)=1;
                    EEG.reject.rejthresh(112)=1;
                    EEG.reject.rejthresh(118)=1;
                    EEG.reject.rejthresh(135)=1;
                    EEG.reject.rejthresh(140:141)=1;
                    EEG.reject.rejthresh(143)=1;
                    EEG.reject.rejthresh(170:171)=1;
                    EEG.reject.rejthresh(186:188)=1;
                    EEG.reject.rejthresh(190)=1;
                    EEG.reject.rejthresh(229)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
                elseif (strcmpi(part_name,'subj27') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(43)=1;
                    EEG.reject.rejthresh(48:49)=1;
                    EEG.reject.rejthresh(51)=1;
                    EEG.reject.rejthresh(61)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);     
                 elseif (strcmpi(part_name,'subj29') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(1)=1;
                    EEG.reject.rejthresh(52)=1;
                    EEG.reject.rejthresh(65)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
                 elseif (strcmpi(part_name,'subj30') && strcmpi(exp.settings,'filt_byTargets_v3'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(113:114)=1;
                    EEG.reject.rejthresh(211:212)=1;
                    EEG.reject.rejthresh(232)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                elseif (strcmpi(part_name,'subj30') && strcmpi(exp.settings,'filt_byCatchTargets_v1_wav'))%only when using this specific settings 
                    EEG.reject.rejthresh = zeros(size(EEG.reject.rejthresh)); %reset variable to all 0s
                    EEG.reject.rejthresh(65)=1;
                    rejtrial(i_set,3).ids = find(EEG.reject.rejthresh==1);
                    EEG = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
                end
                
% ````````````````````````````````````````````````````````````````````````````````````````````

                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % save rejected trials
                EEG.rejtrial = rejtrial;
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                %replace the stored data with this new set
                tempEEG = EEG;               
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                 
                %% Select individual events
                for i_event = 1:nevents
                    EEG = pop_selectevent( tempEEG, 'type', exp.events{i_set,i_event}, 'deleteevents','on','deleteepochs','on','invertepochs','off');
                    EEG = pop_editset(EEG, 'setname', [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}] );
                    EEG = pop_saveset( EEG, 'filename',[part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.set'],'filepath',[exp.pathname '\' exp.setname{i_set} '\Segments\']);
                end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                 
            end %create epoch loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
            %% Time-Frequency Data
            if strcmp('on',exp.tf) == 1 || strcmp('on',exp.singletrials) == 1
                
                if ~exist([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\'])
                    mkdir([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\']);
                end

                
                for i_event = 1:nevents
                    
                    if strcmp('on',exp.epoch) == 0 %loading previous epochs if not created this session
                        filename = [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
                        EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    end
                    
                    if isempty(exp.tfelecs) %if TF electrodes not specified, same as exp.brainelecs
                       exp.tfelecs = exp.brainelecs; 
                    end
                    
                    for i_tf = 1:length(exp.tfelecs)
                        i_chan = exp.tfelecs(i_tf);
                        EEG = eeg_checkset(EEG);
                        [ersp(i_chan,:,:),itc(i_chan,:,:),powbase,times,freqs,dum1,dum2,all_ersp(i_chan).trials] =...
                            pop_newtimef(EEG, 1, i_chan, exp.epochslims*1000, exp.cycles, ...
                            'topovec', i_chan, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo,...
                            'baseline', exp.erspbaseline, 'freqs', exp.freqrange, 'freqscale', 'linear', ...
                            'padratio', exp.padratio, 'plotphase','off','plotitc','off','plotersp','off',...
                            'winsize',exp.winsize,'timesout',exp.timesout);
                    end
                    clear i_chan i_tf

                    if strcmp('on',exp.tf) == 1 %if TF was done already, do not save
                        if exp.cycles(1) > 0 %for wavelet
                            save([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\Wav_' part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'ersp','itc','times','freqs','powbase','exp')
                        else %for FFT
                            save([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\TF_' part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'ersp','itc','times','freqs','powbase','exp')
                        end
                    end
                        
                    
                     % Save single trial data
                    if strcmp('on',exp.singletrials) == 1
                        
                        % Create folder for single trial data
                        if ~exist([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\'],'dir')
                            mkdir([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\']);
                        end
                        
                        % File path name
                        Filepath_Trials = [exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\'];
                        
                        % Save single trial data from the selected electrodes
                        for zzz = 1:length(exp.singletrialselecs)
                            i_chan = exp.singletrialselecs(zzz);
                            elec_all_ersp = all_ersp(i_chan).trials;
                            if exp.cycles(1) > 0 %for wavelet
                                save([Filepath_Trials exp.singtrlelec_name{zzz} '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '_Wav.mat'],...
                                'elec_all_ersp','times','freqs','powbase','exp')
                            else %for FFT
                                save([Filepath_Trials exp.singtrlelec_name{zzz} '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],...
                                'elec_all_ersp','times','freqs','powbase','exp')
                            end
                        end
                        clear i_chan elec_all_ersp zzz
                    end
                    clear Filepath_Trials

                end
                eeglab redraw
                
            end
            clear part_name
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             
        end %i_part loop end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::         
    end %i_set loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    
% !!!!!!!!!!!!!!!!    
catch ME % !!!!!!!
    save('dump') 
    throw(ME)
end % !!!!!!!!!!!!
% !!!!!!!!!!!!!!!! 


% ####################
end % end function ###
% ####################
