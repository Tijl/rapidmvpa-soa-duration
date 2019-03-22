function [] = run_channel_searchlight_multiclass(varargin)

    %%
    if ~ismac %we're on the HPC
        addpath('../CoSMoMVPA/mvpa');
        addpath('../fieldtrip')
    end
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa');
    end
    if isempty(which('ft_defaults'))
        addpath('~/fieldtrip')
    end
    ft_defaults;
    cosmo_warning('off')

    %%
    opt=struct();
    opt=cosmo_structjoin(opt,varargin{:});

    %% 
    if ~isfield(opt,'subject')
        error('subject is required but not specified')
    end
    pp = opt.subject;

    if ~ismac
        % start cluster, give it a unique directory
        % start cluster, give it a unique directory
        % starting a pool can fail when 2 procs are requesting simultaneous
        % thus try again after a second until success
        pool=[];
        while isempty(pool) 
            try
                pc = parcluster('local');
                pc.JobStorageLocation=tempdir;
                pool=parpool(pc);
            catch err
                disp(err)
                delete(gcp('nocreate'));
                pause(1)
            end
        end
        nproc=cosmo_parallel_get_nproc_available();
    end

    %% set and check inputfile
    fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',pp);
    assert(exist(fn,'file')>0,sprintf('file not found: %s',fn))

    %% set outfile
    outfn = sprintf('results/sub-%02i_channel_searchlight_multiclass.mat',pp);

    %% load data
    fprintf('p%i loading\n',pp)
    x=load(fn,'ds');
    cosmo_check_dataset(x.ds,'meeg');
    dsb = cosmo_slice(x.ds,~x.ds.sa.isflipped);
    x=[]; %clear x to save ram
    
    %% make a conditionlist (for easy plotting later)
    conditions=struct();
    for c=1:5
        conditions.number(c) = c;
        x=find(dsb.sa.condition==c,1);
        conditions.durationISI(c) = dsb.sa.durationISI(x);
        conditions.durationSTIM(c) = dsb.sa.durationSTIM(x);
        conditions.durationGAP(c) = dsb.sa.durationGAP(x);
    end
    
    %% make a timevect (for easy plotting later)
    timevect = dsb.a.fdim.values{2};
    
    %% decoding params
    dsb.sa.targets = dsb.sa.stimulusnumber;
    dsb.sa.chunks = dsb.sa.sequencenumber;
    nh1 = cosmo_meeg_chan_neighborhood(dsb, 'count', 4);
    nh2 = cosmo_interval_neighborhood(dsb,'time','radius',0);
    nh = cosmo_cross_neighborhood(dsb,{nh1,nh2});
    m = @cosmo_crossvalidation_measure;
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.progress = 1;
    ma.output='accuracy';
    if ismac
        ma.nproc = 1;
    else
        ma.nproc = 6;
    end
    ma.check_partitions = 0;

    %% decode
    for condition=1:5
        ds = cosmo_slice(dsb,dsb.sa.condition==condition);
        %% set up contrasts
        targets = {ds.sa.levelAnumber,ds.sa.levelBnumber,ds.sa.levelCnumber};
        targetlabels = {'levelA','levelB','levelC'};        
        for t=1:length(targets)
            ds.sa.targets = targets{t};
            fprintf('p%i c%i decoding %s\n',pp,condition,targetlabels{t})
            
            % partitioning scheme
            ma.partitions = cosmo_nfold_partitioner(ds);
            
            % run full searchlight
            r = cosmo_searchlight(ds,nh,m,ma);
                        
            % save under a unique name like res_conditionnumber_levelABC
            eval(sprintf('res_c%i_%s = r;',condition,targetlabels{t}));

            % save all the things
            save(outfn,'res_*','timevect','conditions','-v7.3')
        end
    end
end