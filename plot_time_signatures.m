%% plot_decoding_pairwise
suffix = 'decoding_pairwise_half_sequence';
fprintf('Loading data\n')
files = dir(sprintf('results/boot/stats_%s_*.mat',suffix));
bootstats=struct();
cc = clock();mm='';
for f=1:length(files)
    fn = fullfile(files(f).folder,files(f).name);
    m=load(fn);
    bootstats.MU{f} = m.MU;
    bootstats.BF{f} = m.BF;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end

%% fancy plot including distributions
load('results/stats_decoding_pairwise_half_sequence.mat')

f=figure(1);clf;f.Position=[1 100 1000 600];f.Resize='off';f.PaperPositionMode='auto';
drawnow

% onsets
bfthresh=6;
catnames = {'animacy','category','image'};
onset=[];onset_ci=[];offset=[];offset_ci=[];peak=[];peak_ci=[];
for level = 1:3
    for condi=1:5
        bf = BF{level,condi};
        idx = find(bf>bfthresh);
        x1 = idx(2);
        x2 = idx(end-1);
        
        onset(level,condi) = timevect(x1);
        offset(level,condi) = timevect(x2);
        [~,x] = max(MU{level,condi});
        peak(level,condi) = timevect(x);
        
        boot_x1=[];boot_x2=[];boot_p=[];
        for x = 1:length(bootstats.BF)
            bf = bootstats.BF{x}{level,condi};
            
            idx = find(bf>bfthresh);
            x1 = idx(2);
            x2 = idx(end-1);
            
            boot_x1(x) = x1;
            boot_x2(x) = x2;
            [~,boot_p(x)] = max(bootstats.MU{x}{level,condi});
        end
        onset_ci(level,condi,:) = histcounts(timevect(boot_x1),timevect,'Normalization','pdf');
        offset_ci(level,condi,:) = histcounts(timevect(boot_x2),timevect,'Normalization','pdf');
        peak_ci(level,condi,:) = histcounts(timevect(boot_p),timevect,'Normalization','pdf');
        
    end
end

co = viridis(6);
co2 = plasma(4);
aw=.05;
mw = 5;
mscale = 15;
for level=1:3
    a=subplot(2,3,level);hold on
    title(sprintf('decoding window - %s',catnames{level}))
    lb=[];
    for condi=1:5
        y = 6-condi;
        x = squeeze(onset_ci(level,condi,:));
        fill([timevect fliplr(timevect)],y+aw+mscale*[movmean(x,mw)' 0 0*timevect],co(condi,:),'FaceAlpha',.2,'EdgeAlpha',0);hold on
        x = squeeze(offset_ci(level,condi,:));
        fill([timevect fliplr(timevect)],y-aw-mscale*[movmean(x,mw)' 0 0*timevect],co(condi,:),'FaceAlpha',.2,'EdgeAlpha',0);hold on
        x = [onset(level,condi) offset(level,condi)];
        fill([x fliplr(x)],y+[-aw -aw aw aw],co(condi,:),'EdgeAlpha',0);hold on
        lb{condi} = sprintf('duration: %i SOA: %i',1000*conditions.durationSTIM(condi),1000*conditions.durationISI(condi));
        text(x(1),y,sprintf('%ims',x(1)),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',9,'Color',co(condi,:))
        text(x(2),y,sprintf('%ims',x(2)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9,'Color',co(condi,:))
    end
    xlabel('time')
    a.YTick=1:5;
    a.YLim=[0.5 5.5];
    a.XLim=[-100 1000];
    a.XTick=0:200:1000;
    a.YTickLabel=fliplr(lb);
end


%peak time
co = viridis(6);
co2 = plasma(4);
aw=.2;
for level=1:3
    a=subplot(2,3,3+level);hold on
    title(sprintf('peak latency - %s',catnames{level}))
    lb=[];
    for condi=1:5
        y = 6-condi;
        x = squeeze(peak_ci(level,condi,:))';
        fill([timevect fliplr(timevect)],y+mscale*[movmean(x,mw) 0 0*timevect],co(condi,:),'FaceAlpha',.2,'EdgeAlpha',0);hold on
        x = peak(level,condi);
        plot(x,y,'.','Color',co(condi,:),'MarkerSize',10);hold on
        lb{condi} = sprintf('duration: %i SOA: %i',1000*conditions.durationSTIM(condi),1000*conditions.durationISI(condi));
        text(x,y,sprintf('%ims',x),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',9,'Color',co(condi,:))
    end
    xlabel('time (ms)')
    a.YTick=1:5;
    a.XLim=[-100 1000];
    a.XTick=0:200:1000;
    a.YTickLabel=fliplr(lb);
    a.YLim=[0.5 5.5];
end

fn = 'figures/figure_time_signatures_distributions';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
