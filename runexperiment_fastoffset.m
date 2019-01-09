%% Fast presentation experiment with varying durations and ISI

%% set up psychtoolbox if needed (should only be needed for testing)
if isempty(which('PsychtoolboxVersion.m'))
    addpath(genpath('~/PHD/Matlab/Psychtoolbox/'))
end
addpath(genpath('SlackMatlab'));

close all; %close all figure windows
clearvars; %clear all variables

%% START define parameters
p=struct();

%testing parameters: CHECK THESE!
p.isrealexperiment = 1; %should only be 0 for testing; skips overwrite checks and doesn't ask for subject nr
p.fullscreen = 1; %should only be 0 for testing; if 0 it does not create a fullscreen window
p.isEEGexperiment = 1; %should only be 0 for testing; does not send actual triggers
p.synctest = 1; %should only be 0 for testing; skips synchronisation tests

% trigger port:
p.triggerportcode = 'CFF8';
p.triggerstreamstart = 1;
p.triggerstimon = 2;
p.triggerstimoff = 3;
p.triggerquerystimon = 4;
p.triggerquerystimoff = 5;
p.triggerduration = .01;

% stimulus parameters
p.nstimuli = 24; % number of image stimuli, should be in ./Stimuli
p.nsequencespercondition = 40; %number of sequences per condition (40 in real)
% p.nsequencespercondition = 10; %number of sequences per condition (10 for students)
p.nconditions = 5;
p.nsequences = p.nsequencespercondition*p.nconditions;
p.durationISI  = [.2 .2 .2  .1  .05]; %isi for each condition
p.durationSTIM = [.2 .1 .05 .05 .05]; %duration for stimulus
p.durationGAP = p.durationISI-p.durationSTIM; %duration for gap

%other general timing parameters
p.sequencefixduration = 1; %time before first stimulus in a sequence;
p.sequencepostfixduration = 1; %time after last stimulus in a sequence;
p.querystimduration = .2; %duration of query stimulus
p.feedbackduration = .5; %duration of feedback to query

%display parameters
p.stimulussize=[128 128]; %stimulus size onscreen (pixels)
p.fixationsize=20; %diameter of fixation cross (pixels)
p.fixationcolour=[0 0 0]; %colour of fixation cross
p.backgroundcolor=[127 127 127]; %gray
p.fontsize=22;
p.windowsize=[0 0 800 600]; %only used if not running fullscreen

%slack channel
p.strHookURL = 'https://hooks.slack.com/services/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';

%% END define parameters


%% get subject info
if ~p.isrealexperiment
    p.subjectnr = 0;
    disp('    +-----------------------------------------+')
    disp('    | WARNING: p.isrealexperiment is set to 0 |')
    disp('    | If this is a real EEG run, set to 1     |')
    disp('    +-----------------------------------------+')
    WaitSecs(1);
else
    p.subjectnr = input('\n    Subject number: ','s');
    p.subjectnr = str2double(p.subjectnr);
    if isnan(p.subjectnr)
        error('ERROR: invalid subject number');
    end
end
p.datafilename = sprintf('data_subject_%i.mat',p.subjectnr);

%check if we are possibly overwriting data
if p.isrealexperiment
    if exist(p.datafilename,'file')
        error(['ERROR: ' p.datafilename ' exists. Move (or delete) this file or use a different subject number']);
    end
end

%This should only be used to test/debug. CHECK THIS!!
if p.synctest
    Screen('Preference', 'SkipSyncTests', 0);
else
    Screen('Preference', 'SkipSyncTests', 1);
    disp('    +-------------------------------------------------+')
    disp('    | WARNING: p.synctest is set to 0                 |')
    disp('    | If this is an EEG experiment, p.synctest to 1   |')
    disp('    +-------------------------------------------------+')
    WaitSecs(1);
end

% set up seed
p.randomseed = rng(p.subjectnr);

p.responsekeyswap = ~mod(p.subjectnr,2);
if p.responsekeyswap
    p.leftresponsename = 'Animal';
    p.rightresponsename = 'Vehicle';
else
    p.leftresponsename = 'Vehicle';
    p.rightresponsename = 'Animal';
end

%% open i/o port
if p.isEEGexperiment
    %create io32 interface object
    p.ioObj = io64;
    % check the port status
    status = io64(p.ioObj);
else
    disp('    +-------------------------------------------------+')
    disp('    | WARNING: p.isEEGexperiment is set to 0          |')
    disp('    | If this is an EEG experiment, set to 1          |')
    disp('    +-------------------------------------------------+')
    WaitSecs(1);
    status = 1;
end
if status == 0
    p.address = hex2dec(p.triggerportcode); %stim2
    display(['Functioning parallel port opened at: ' num2str(p.address)])
else
    p.ioObj = [];
    p.address = [];
end

%make a copy of the raw experiment code and store it so we can recreate the exact script
p.experimentcode = fileread('runexperiment_fastoffset.m');

%% load stimuli

p.stimuli = struct('raw',[],'alpha',[],'rawalpha',[],'name',[]);

for i=1:p.nstimuli
    
    imfiles=dir(sprintf('Stimuli_resized/stim%03i_*.png',i));
    if length(imfiles)>1
        error(sprintf('Multiple .png files start with %i_ in Stimuli',i)) %#ok<SPERR>
    end
    
    imname = imfiles(1).name;
    impath = ['Stimuli_resized' filesep imname];
    [im, ~, alpha]=imread(impath);
    
    parts=split(imname,{'_','.'}); 
    
    p.stimuli(i).raw = double(im);
    p.stimuli(i).alpha = double(alpha);
    p.stimuli(i).rawalpha = cat(3,p.stimuli(i).raw,p.stimuli(i).alpha);
    p.stimuli(i).path = impath;
    p.stimuli(i).name = imname;
    p.stimuli(i).animal = strcmpi(parts{3},'animals');
    p.stimuli(i).stimulusnumber = i;
end

%% create trialstructure
[~,X]=sort(rand(p.nstimuli,p.nsequences));
%padding
p.ispadding = [true(p.nstimuli/2,p.nsequences); false(p.nstimuli,p.nsequences); true(p.nstimuli/2,p.nsequences)];
p.trialstructure = zeros(size(p.ispadding));
p.trialstructure(p.ispadding) = flipud(X);
p.trialstructure(~p.ispadding) = X;
% create conditionstructure
[~,X] = sort(rand(p.nconditions,p.nsequencespercondition));
p.conditionlist = X(:);
% create querystimuli
[~,X] = sort(rand(p.nstimuli,p.nsequences));
p.querystimuli = fliplr(X(1:p.nsequences)); %flip so can remove the first few
p.querystimulianimal = [p.stimuli(p.querystimuli).animal];

%% setup abort function
abort = @(x) eval(['sca;Priority(0);ListenChar(0);ShowCursor();fprintf(''\n\nABORTING... emergency save...\n'');save(''aborted.mat'');sca;'...
    'error(struct(''stack'',struct(''file'',''runexperiment_fastoffset'',''name'',''runexperiment_fastoffset'',''line'',0),'...
    '''message'',''#### EXPERIMENT ABORTED BY USER ####''));']);

%% disable keyboard input
ListenChar(2);
KbName('UnifyKeyNames') %cause old version of psychtoolbox
KbCheck(); % make sure kbcheck is compiled so it is always accurate on first call
if p.responsekeyswap
    p.key_animal = KbName('s');
    p.key_vehicle = KbName('d');
else
    p.key_animal = KbName('d');
    p.key_vehicle = KbName('s');
end

%% summarize parameters
disp('    +------------------------+')
disp('    | Experiment parameters: |')
disp('    +------------------------+')
WaitSecs(1);
disp(p);
disp('If this is correct, press return to start.')
[~, keycode, ~] = KbWait(); %wait for a keypress
if ~keycode(KbName('Return'));abort();end % check for return key

%% open window, and wait for the photodiode setup
PsychDefaultSetup(2);
p.screennumber=max(Screen('Screens'));

p.black=BlackIndex(p.screennumber);
p.gray=GrayIndex(p.screennumber);
p.white=WhiteIndex(p.screennumber);
if p.fullscreen
    p.windowsize=[];
elseif ismac && max(Screen('Screens',1))==2
    p.windowsize=[];
end
[p.window, p.windowrect]=Screen('OpenWindow',p.screennumber, p.backgroundcolor,p.windowsize,32,2);
Screen('TextSize',p.window,p.fontsize);
Screen('BlendFunction',p.window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

p.flipinterval = Screen('GetFlipInterval',p.window);
p.flipadjust = .5*p.flipinterval;

%fixation fuction
p.fixationlocation = .5*p.windowrect([3 3 3 3;4 4 4 4])+.5*[-p.fixationsize p.fixationsize 0 0;0 0 -p.fixationsize p.fixationsize];
drawfixation = @() Screen('DrawLines', p.window, p.fixationlocation,2,p.fixationcolour);
drawfixationred = @() Screen('DrawLines', p.window, p.fixationlocation,2,[200 0 0]);
drawfixationgreen = @() Screen('DrawLines', p.window, p.fixationlocation,2,[0 200 0]);

%% make textures
textures = zeros(1,p.nstimuli);
texturesflipped = zeros(1,p.nstimuli);
for stim = 1:p.nstimuli
    textures(stim) = Screen('MakeTexture', p.window, p.stimuli(stim).rawalpha);
    texturesflipped(stim) = Screen('MakeTexture', p.window, fliplr(p.stimuli(stim).rawalpha));
end

%% start sequence loop
data = struct();
data.imagename = {};
[data.time_stimon, data.time_stimoff, ...
data.isflipped, data.stimulusnumber, ...
data.condition, data.durationGAP, ...
data.durationSTIM, data.durationISI, ...
data.sequencenumber, data.presentationnumber, ...
data.time_sequencestart, data.querystim, ...
data.time_querystimon, data.time_querystimoff, ...
data.response, data.correct, data.rt] = deal([]);

for sequencenumber = 1:p.nsequences
    condition = p.conditionlist(sequencenumber);
    durationSTIM = p.durationSTIM(condition);
    durationGAP = p.durationGAP(condition);
    durationISI = p.durationISI(condition);
    
    DrawFormattedText(p.window, sprintf('%i / %i\n\nPress any button',sequencenumber,p.nsequences),'center', 'center', p.white);
    Screen('Flip', p.window);
    %wait for keypress to continue
    [~, keycode] = KbWait();
    if keycode(KbName('escape'));abort();end
    
    drawfixation();
    time_fixon = Screen('Flip', p.window);
    data.time_sequencestart(end+1)=time_fixon;
    if p.isEEGexperiment
        io64(p.ioObj, p.address, p.triggerstreamstart);
        WaitSecs(p.triggerduration);
        io64(p.ioObj, p.address, 0);
    end
    drawfixation();
    for presentationnumber = 1:size(p.trialstructure,1)
        stim = p.trialstructure(presentationnumber,sequencenumber);
        ispadding = p.ispadding(presentationnumber,sequencenumber);
        
        data.sequencenumber(end+1) = sequencenumber;
        data.presentationnumber(end+1) = presentationnumber;
        data.imagename{end+1}=p.stimuli(stim).name;
        data.stimulusnumber(end+1)=stim;
        data.isflipped(end+1)=ispadding;
        data.durationGAP(end+1) = durationGAP;
        data.durationSTIM(end+1) = durationSTIM;
        data.durationISI(end+1) = durationISI;
        data.condition(end+1) = condition;
        
        if ispadding
            tex = texturesflipped(stim);
        else
            tex = textures(stim);
        end
        Screen('DrawTexture',p.window,tex,[],CenterRect([0 0 p.stimulussize],p.windowrect));
        drawfixation();
        if presentationnumber==1
            data.time_stimon(end+1) = Screen('Flip', p.window,time_fixon+p.sequencefixduration-p.flipadjust);
        else
            data.time_stimon(end+1) = Screen('Flip', p.window, data.time_stimon(end)+durationGAP+durationSTIM-p.flipadjust);
        end
        if p.isEEGexperiment
            io64(p.ioObj, p.address, p.triggerstimon);
            WaitSecs(p.triggerduration);
            io64(p.ioObj, p.address, 0);
        end
        if durationGAP>0
            drawfixation();
            data.time_stimoff(end+1) = Screen('Flip', p.window, data.time_stimon(end)+durationSTIM-p.flipadjust);        
            if p.isEEGexperiment
                io64(p.ioObj, p.address, p.triggerstimoff);
                WaitSecs(p.triggerduration);
                io64(p.ioObj, p.address, 0);
            end
        elseif presentationnumber>1
            data.time_stimoff(end+1) = data.time_stimon(end);
        end
    end
    if durationGAP==0
        drawfixation();
        data.time_stimoff(end+1) = Screen('Flip', p.window, data.time_stimon(end)+durationSTIM-p.flipadjust);
    end
    WaitSecs(p.sequencepostfixduration);
    
    %now do the query + response
    stim = p.querystimuli(sequencenumber);
    tex = textures(stim);
    Screen('DrawTexture',p.window,tex,[],CenterRect([0 0 p.stimulussize],p.windowrect));
    drawfixation();
    
    data.time_querystimon(end+1) = Screen('Flip', p.window);
    if p.isEEGexperiment
        io64(p.ioObj, p.address, p.triggerquerystimon);
        WaitSecs(p.triggerduration);
        io64(p.ioObj, p.address, 0);
    end

    rect = CenterRect([0 0 p.stimulussize],p.windowrect);    
    DrawFormattedText(p.window, p.leftresponsename, 'right', 'center', p.white, [], [], [], [], [], [p.windowrect(1:2) rect(1)-20 p.windowrect(4)]);
    DrawFormattedText(p.window, p.rightresponsename, rect(3)+20, 'center', p.white, [], [], [], [], []);
    
    drawfixation();
    data.time_querystimoff(end+1) = Screen('Flip', p.window, data.time_querystimon(end)+p.querystimduration-p.flipadjust);
    if p.isEEGexperiment
        io64(p.ioObj, p.address, p.triggerquerystimoff);
        WaitSecs(p.triggerduration);
        io64(p.ioObj, p.address, 0);
    end
    
    response = -1; keypressed=0;
    % simulate responses
    % response = rand>.5; rt = .5*rand+.4; WaitSecs(rt);
    
    while response<0
        [keypressed, secs, keycode] = KbCheck();
        if keycode(KbName('escape'));abort();end
        keys = keycode([p.key_vehicle p.key_animal]);
        if any(keys)
            response = find(keys)-1;
            rt = secs-data.time_querystimon(end);
        end
    end
    correct = response == p.querystimulianimal(sequencenumber);
    data.querystim(end+1) = stim;
    data.response(end+1) = response;
    data.rt(end+1) = rt;
    data.correct(end+1) = correct;
    if correct
        drawfixationgreen();
    else
        drawfixationred();
    end
    Screen('Flip', p.window);
    WaitSecs(p.feedbackduration);
    drawfixation();
    Screen('Flip', p.window);
    % save data every n sequences
    if ~mod(sequencenumber,p.nconditions)
        save(p.datafilename,'p','data');
        
        try
            %slack bot
            message = sprintf('<@tijlgrootswagers> and <@amanda>, pp%i seq %i/200 done',p.subjectnr,sequencenumber);
            SendSlackNotification(p.strHookURL,message,'#eeglab','MATLAB');
        catch
        end
    end
    % wait for key release
    while keypressed
        keypressed = KbCheck();
    end
end

%% save
save(p.datafilename,'p','data');

%% close all
Priority(0);ListenChar(0);ShowCursor();
Screen('CloseAll');

%% plot timing data
% figure(1);plot([[0 diff(data.time_stimon)]; data.time_stimoff-data.time_stimon]');ylim([0 .5]);legend({'ISI','DUR'})
