function [] = run_decoding_pairwise_half_sequence(varargin)

    %%
    if ~ismac %we're on the HPC
        addpath('../CoSMoMVPA/mvpa');
    end
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa');
    end

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
    outfn = sprintf('results/sub-%02i_decoding_pairwise_half_sequence.mat',pp);

    %% load data
    fprintf('p%i loading\n',pp)
    x=load(fn,'ds');
    cosmo_check_dataset(x.ds,'meeg');
    dsb = cosmo_slice(x.ds,~x.ds.sa.isflipped);
    x=[]; %clear data to save ram
    
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
    nh = cosmo_interval_neighborhood(dsb,'time','radius',0);
    m = @cosmo_crossvalidation_measure;
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.progress = 1;
    ma.output='fold_accuracy';
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
        cvtargets = {ds.sa.levelCnumber,ds.sa.levelCnumber,ds.sa.levelCnumber,ds.sa.levelCnumber};
        targetlabels = {'levelA','levelB','levelC'};        
        for t=1:length(targets)
            ds.sa.targets = targets{t};
            ds.sa.cvtargets = cvtargets{t};
            fprintf('p%i c%i decoding %s\n',pp,condition,targetlabels{t})
            % all pairwise combinations
            combs = combnk(unique(ds.sa.targets,'rows'),2);
            % all chunks to leave out
            uc = unique(ds.sa.chunks);
            % create exemplar-by-sequence partitioning scheme
            ma.partitions = struct();
            ma.partitions.train_indices = {};
            ma.partitions.test_indices = {};
            sa=struct('target1',[],'target2',[],'leftoutchunk',[],'leftoutexemplar1',[],'leftoutexemplar2',[]);
            for i=1:size(combs,1) % for each pairwise comparision to test
                % find the epochs in this pair
                idx1 = ismember(ds.sa.targets,combs(i,:));
                % ue1 and ue2 are the unique exemplars (to leave out in the test set)
                ue1 = unique(ds.sa.cvtargets(ds.sa.targets==combs(i,1)));
                ue2 = unique(ds.sa.cvtargets(ds.sa.targets==combs(i,2)));
                % leave all combinations of exemplar pairs out once
                ue = [repelem(ue1,length(ue2),1) repmat(ue2,length(ue1),1)];
                for j=1:length(uc) % for each chunk to leave out
                    idx2 = ds.sa.chunks==uc(j); % find chunk to leave out
                    for k=1:size(ue,1) % for each exemplar pair to leave out
                        % store targets in results
                        sa.target1(end+1,1) = combs(i,1);
                        sa.target2(end+1,1) = combs(i,2);
                        % store left out chunk and exemplar in result
                        sa.leftoutchunk(end+1,1) = uc(j);
                        sa.leftoutexemplar1(end+1,1) = ue(k,1);
                        sa.leftoutexemplar2(end+1,1) = ue(k,2);
                        % set partitions
                        % if size(ue,1)>1 then we are doing
                        % exemplar-by-sequence (otherwise just sequence,
                        % (for the lowest (image) level)
                        if size(ue,1)>1
                            idx3 = ismember(ds.sa.cvtargets,ue(k,:));
                            ma.partitions.train_indices{1,end+1} = find(idx1 & ~idx2 & ~idx3);
                            ma.partitions.test_indices{1,end+1} = find(idx1 & idx2 & idx3);
                        else
                            ma.partitions.train_indices{1,end+1} = find(idx1 & ~idx2);
                            ma.partitions.test_indices{1,end+1} = find(idx1 & idx2);
                        end
                    end
                end
            end
            % run full searchlight
            r = cosmo_searchlight(ds,nh,m,ma);
            % merge fold information into result (targets, left out chunk & exemplar)
            r.sa = cosmo_structjoin(r.sa,sa);
            
            % save under a unique name like res_conditionnumber_levelABC
            eval(sprintf('res_c%i_%s = r;',condition,targetlabels{t}));

            % save all the things
            save(outfn,'res_*','timevect','conditions','-v7.3')
        end
    end
end