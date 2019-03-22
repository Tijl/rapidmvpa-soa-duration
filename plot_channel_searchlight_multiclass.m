
%%
res_cell={};
for i = 1:20
    x = load(sprintf('results/sub-%02i_channel_searchlight_multiclass.mat',i));
    for c = 1:5
        for t=1:3
            res = x.(sprintf('res_c%i_level%s',c,char('A'+t-1)));
            res.sa.c = c;
            res.sa.t = t;
            res_cell{end+1} = res;
        end
    end
    conditions = x.conditions;
    timevect = x.timevect;
end
res_all = cosmo_stack(res_cell);


%%
layout=cosmo_meeg_find_layout(res);

plotnr=0;
f=figure(1);clf
drawnow
f.Position=[f.Position(1:2) 1200 1600];
f.Resize='off';
f.PaperPositionMode='auto';f.PaperOrientation='portrait';
levelnames = {'animacy','object','image'};

h0mean = [1/2 1/6 1/24];
for c=1:5
    for t=1:3
        plotnr=plotnr+1;
        subplot(5,3,plotnr)
        
        x = cosmo_slice(res_all,res_all.sa.c==c & res_all.sa.t==t);
        x.samples = x.samples-h0mean(t);
        
        % map to FT struct for visualization
        ft=cosmo_map2meeg(x);

        co = viridis(5);
        % show figure with plots for each sensor
        cfg = [];
        cfg.interactive = 'yes';
        cfg.ylim        = 'maxmin';
        cfg.layout      = layout;
        cfg.showcomment = 'no';
        cfg.graphcolor  = co(1,:);
        cfg.showscale = 'no';
        ft_multiplotER(cfg, ft);
        
        tt = sprintf('duration: %i soa: %i - %s', 1000*conditions.durationSTIM(c), 1000*conditions.durationISI(c),levelnames{t});
        
        title(tt)
        drawnow
    end
end

%% save
fn = 'figures/figure_channel_searchlight_multiclass';
f.PaperOrientation='portrait';
print(f,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');


%% 2 time points (peaks) and plot per level
layout=cosmo_meeg_find_layout(res);

h0mean = [1/2 1/6 1/24];
ranges = [.025 .04 .015];
levelnames = {'animacy','object','image'};
timewindows = {[100 150],[175 225]};
for t=1:3
    
    plotnr=0;
    f=figure(t);clf
    drawnow
    f.Position=[f.Position(1:2) 500 1000];
    f.Resize='off';
    f.PaperPositionMode='auto';f.PaperOrientation='portrait';

    for c=1:5
        for tw=1:2
            plotnr=plotnr+1;
            subplot(5,2,plotnr)

            x = cosmo_slice(res_all,res_all.sa.c==c & res_all.sa.t==t);
            x.samples = x.samples-h0mean(t);

            % map to FT struct for visualization
            ft=cosmo_map2meeg(x);
            ft = ft_timelockanalysis([],ft);

            co = viridis(5);
            % show figure with plots for each sensor
            cfg = [];
            cfg.zlim        = [0 ranges(t)];
            cfg.layout      = layout;
            cfg.comment = 'no';
            %cfg.graphcolor  = co(1,:);
            cfg.showscale = 'no';
            cfg.xlim = timewindows{tw};
            cfg.style = 'straight';
            cfg.colormap = viridis();
            ft_topoplotER(cfg, ft);
            
            drawnow
            
            a = gca;p=a.Position;
            cc = colorbar;
            a.Position = p;
            
            cc.Ticks=[];
            if tw==2
                cc.Visible='off';
            else
                cc.Ticks = [0 1].*ranges(t);
                cc.TickLabels = {sprintf('%.2f (chance)',100*h0mean(t)),sprintf('%.2f',100*(h0mean(t)+ranges(t)))};
            end
                
            
            tt = {sprintf('%s - [%i:%i]ms', levelnames{t},timewindows{tw}),...
                sprintf('duration: %i soa: %i', 1000*conditions.durationSTIM(c), 1000*conditions.durationISI(c))};

            %title(tt,'Units','Normalized','Position',[0,1,0],'HorizontalAlignment','left','VerticalAlignment','bottom')
            title(tt)
            drawnow
        end
    end
    
    fn = sprintf('figures/figure_channel_searchlight_multiclass_%s',levelnames{t});
    f.PaperOrientation='portrait';
    print(f,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
end

