function stats_decoding_pairwise_half_sequence()
    
    %%
    if ~ismac %we're on the HPC
        addpath('../CoSMoMVPA/mvpa');
    end
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa');
    end

    %% load data
    suffix = 'decoding_pairwise_half_sequence';
    fprintf('Loading data\n')
    nsubjects = length(dir(sprintf('results/sub-*%s.mat',suffix)));
    res_cell={};
    cc = clock();mm='';
    for f=1:nsubjects
        fn = sprintf('results/sub-%02i_%s.mat',f,suffix);
        m=load(fn);
        tn = {'A','B','C'};
        for t=1:3
            for c=1:5
                try
                    eval(sprintf("x = cosmo_average_samples(m.res_c%i_level%s,'split_by',{});",c,tn{t}));
                    x.sa={};
                    x.sa.level = t*ones(size(x.samples,1),1);
                    x.sa.levelname = repmat(char(t+'A'-1),size(x.samples,1),1);
                    x.sa.subject = f*ones(size(x.samples,1),1);
                    x.sa.condition = c*ones(size(x.samples,1),1);
                    res_cell{t,c,f} = x;
                catch
                    fprintf('MISSING: m.res_c%i_level%s\n',c,tn{t})
                end
            end
        end
        timevect = m.timevect;
        conditions = m.conditions;
        mm = cosmo_show_progress(cc,f/nsubjects,sprintf('%i/%i',f,nsubjects),mm);
    end
    
    %%
    rng(1); %replicable bootstrap sampling
    mmboot='';ccboot=clock();
    nboot=1000;
    for bootnr = 0:nboot
    
        if bootnr
            bootpp = randi(20,1,20);
            res_all = cosmo_stack(res_cell(:,:,bootpp));
        else
            fprintf('bootnr %i\n',bootnr)
            res_all = cosmo_stack(res_cell);
        end
        
        %% BF against chance
        if ~bootnr
            fprintf('BF against chance\n')
        end
        combs=[];
        for t=1:3
            for c=1:5
                combs(end+1,1:2) = [t,c];
            end
        end
        [MU,SE,TSTAT,BF,BFT] = deal(cell(3,5));
        cc = clock();mm='';
        for k=1:length(combs) 
            t = combs(k,1);
            c = combs(k,2);
            idx = res_all.sa.condition==c & res_all.sa.level==t;
            x = res_all.samples(idx,:);
            if isempty(x)
                x = zeros(2,size(res_all.samples,2));
            elseif size(x,1)==1
                x = [x;x];
            end
            mu = mean(x);
            se = std(x)./sqrt(nsubjects);
            MU{t,c} = mu;
            SE{t,c} = se;
            ts = (mu-.5)./se;
            TSTAT{t,c} = ts;
            tidx = timevect<0;
            tlower = max(mu(tidx)-.5);
            bf = zeros(1,length(timevect));
            for i=1:length(timevect)
                bf(i) = bayesfactor(mu(i)-.5,se(i),1,tlower,.5);
            end
            bft = t1smpbf(ts,nsubjects);
            BF{t,c} = bf;
            BFT{t,c} = bft;
            %mm = cosmo_show_progress(cc,k/length(combs),sprintf('%i/%i (%i,%i) ',k,length(combs),t,c),mm);
        end

        %% BF difference
        if ~bootnr
            fprintf('BF difference\n')
            [MUdiff,SEdiff,BFdiff,TSTATdiff,BFTdiff] = deal({});
            combs=[];
            for t=1:3
                for c1=1:5
                    for c2=1:5
                        if c2<=c1;continue;end
                        combs(end+1,1:3)=[t,c1,c2];
                    end
                end
            end
            cc = clock();mm='';
            for k = 1:length(combs)
                t = combs(k,1);
                c1 = combs(k,2);
                c2 = combs(k,3);
                idx1 = res_all.sa.condition==c1 & res_all.sa.level==t;
                x = res_all.samples(idx1,:);
                idx2 = res_all.sa.condition==c2 & res_all.sa.level==t;
                y = res_all.samples(idx2,:);
                mu = mean(x-y);
                se = std(x-y)./sqrt(nsubjects);
                MUdiff{t,c1,c2} = mu;
                MUdiff{t,c2,c1} = mu;
                SEdiff{t,c1,c2} = se;
                SEdiff{t,c2,c1} = se;
                ts = mu./se;
                TSTATdiff{t,c1,c2} = ts;
                TSTATdiff{t,c2,c1} = ts;
                tidx = timevect<0;
                tlower = max(abs(mu(tidx)));
                bf = zeros(1,length(timevect));
                for i=1:length(timevect)
                    bf(i) = bayesfactor(abs(mu(i)),se(i),1,tlower,.5);
                end
                bft = t1smpbf(ts,nsubjects);
                BFdiff{t,c1,c2} = bf;
                BFdiff{t,c2,c1} = bf;
                BFTdiff{t,c1,c2} = bft;
                BFTdiff{t,c2,c1} = bft;
                mm = cosmo_show_progress(cc,k/length(combs),sprintf('%i/%i (%i,%i) ',k,length(combs),t,c),mm);
            end
        end
        if bootnr
            save(sprintf('results/boot/stats_decoding_pairwise_half_sequence_boot%i.mat',bootnr),'MU*','SE*','BF*','TSTAT*','timevect','conditions','-v7.3')
            mmboot = cosmo_show_progress(ccboot,bootnr/nboot,sprintf('%i/%i',bootnr,nboot),mmboot);
        else
            save('results/stats_decoding_pairwise_half_sequence.mat','MU*','SE*','BF*','TSTAT*','timevect','conditions','-v7.3')
        end
    end