function preprocessing_fastoffset(pp)

%% cosmo
if ~ismac
    addpath('../CoSMoMVPA/mvpa')
end

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% eeglab
if ~ismac
    addpath('../eeglab')
end

if isempty(which('eeglab'))
    addpath('~/Dropbox (Personal)/MATLAB/eeglab14_0_0b')
end
eeglab;

%% load eventlist
eventfn = sprintf('data/sub-%02i/eeg/sub-%02i_task-rsvp_events.tsv',pp,pp);
tmp = tempname;
copyfile(eventfn,tmp)
eventlist=readtable(tmp);

%% check for preprocessed continuous already
contfn = sprintf('data/derivatives/eeglab/sub-%02i_task-rsvp_continuous.set',pp);
if exist(contfn,'file')
    EEG_cont=pop_loadset(contfn)
else
    % load EEG file
    EEG_raw = pop_loadbv(sprintf('data/sub-%02i/eeg/',pp),sprintf('sub-%02i_task-rsvp_eeg.vhdr',pp));
    EEG_raw = eeg_checkset(EEG_raw);
    EEG_raw.setname = sprintf('sub-%02i',pp);
    EEG_raw = eeg_checkset(EEG_raw);

    % high pass filter
    EEG_raw = pop_eegfiltnew(EEG_raw, 0.1,[]);
    
    % low pass filter
    EEG_raw = pop_eegfiltnew(EEG_raw, [],100);
    
    % downsample
    EEG_raw = pop_resample(EEG_raw, 250);
    EEG_raw = eeg_checkset(EEG_raw);
    
    % remove extra triggers for participants with script crash
    if pp==8
        % get trigger latencies and find big break
        lat = [EEG_raw.event(:).latency];
        part2start = find(diff(lat)==max(diff(lat)))+1; % trigger num of second part (after break)
        % change previous sequence's triggers to be 9s instead
        nbad = 48+1+2; % for each trial, 48 images plus fixation trigger plus two trigs at end
        badstart = part2start-nbad;
        for n = badstart:(part2start-1)
            EEG_raw.event(n).type = 9;
        end
        EEG_raw = eeg_checkset(EEG_raw);
    elseif pp==18
        % keep sequences 1-170 and last 30
        trigs = {EEG_raw.event(:).type};
        crittrigs = find(strcmp(trigs,'E  1')); % start of each sequence
        
        % mark triggers signifying not in 1-170 or last 30 as '9'
        badtrigs = crittrigs(171):crittrigs(length(crittrigs)-29);
        for n = 1:length(badtrigs)
            EEG_raw.event(badtrigs(n)).type = 9;
        end
        EEG_raw = eeg_checkset(EEG_raw);
    end
    
    % create eventlist and save
    EEG_cont = pop_creabasiceventlist( EEG_raw , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' });
    EEG_cont = eeg_checkset(EEG_cont);
    
    pop_saveset(EEG_cont,contfn);
end
    
% run binlister for bins and extract T2 bin-based epochs
EEG_cont = pop_binlister( EEG_cont , 'BDF', 'fastoffsets_makebins.txt', 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'ExportEL', tempname);
EEG_cont = pop_epochbin( EEG_cont , [-100  1000], 'none'); % no baseline correct
EEG_cont = eeg_checkset(EEG_cont);

%% convert to cosmo
ds = cosmo_flatten(permute(EEG_cont.data,[3 1 2]),{'chan','time'},{{EEG_cont.chanlocs.labels},EEG_cont.times},2);
ds.a.meeg=struct(); %or cosmo thinks it's not a meeg ds 
ds.sa = table2struct(eventlist,'ToScalar',true);
cosmo_check_dataset(ds,'meeg');

% save cosmo ds
save(sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',pp), 'ds', '-v7.3')
    
