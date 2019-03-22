%% plot_decoding_pairwise
load('results/stats_decoding_pairwise_half_sequence.mat')

f=figure(1);clf;f.Position=[1 100 1000 1000];f.Resize='off';f.PaperPositionMode='auto';
drawnow
co = viridis(6);
co2 = plasma(4);
bfthresh = 6;
st = .48; %upper limit
lh = .003; %spacing between dot lines
li = -.012; %spacing between triplets
ms = 2; %dot size
zc = .9*[1 1 1]; %empty space bf colour
catnames = {'Animacy','Object','Image'};
for level = 1:3
    
    bn = {'200ms','100ms','50ms'};
    bnd = {'200ms vs 100ms','200ms vs 50ms','100ms vs 50ms'};
    bd = {BFdiff{level,1,2} BFdiff{level,1,3} BFdiff{level,2,3};...
          BFdiff{level,3,4} BFdiff{level,3,5} BFdiff{level,4,5}};
    ttext = {sprintf('%s - Effect of duration (SOA: 200ms)',catnames{level}),...
        sprintf('%s - Effect of SOA (duration: 50ms)',catnames{level})};
    dn = {conditions.durationSTIM,conditions.durationISI};
    dt = {'duration','SOA'};
    idx = [1:3;3:5];
    
    for plotnr = 1:2  
        a=subplot(3,2,2*level+plotnr-2);hold on  
        plot(timevect,.5+0*timevect,'k--')
        title(ttext{plotnr})
        h=[];
        
        if plotnr==1
            %connecting lines
            for z=-1:2
                i=3;
                plot([-190 min(timevect)],[.435+z*.015 (st+lh*z+li*(i-1))],'-o','Color',zc.*.9,'LineWidth',ms,'Clipping','off','MarkerSize',ms,'MarkerFaceColor',zc*.9);
            end
        end
        
        for i=1:3
            ii=idx(plotnr,i);
            mu = MU{level,ii};
            se = SE{level,ii};
            bf = (BF{level,ii});
            fill([timevect,fliplr(timevect)],[mu-se fliplr(mu+se)],co(ii,:),'FaceAlpha',.2,'EdgeAlpha',0)
            h(i)=plot(timevect,mu,'-','Color',co(ii,:),'LineWidth',1,...
                'DisplayName',[sprintf('%s: %ims',dt{3-plotnr},1000*dn{3-plotnr}(ii)),...
                ' & ' sprintf('%s: %ims',dt{plotnr},1000*dn{plotnr}(ii))]);

            x = zeros(size(bf));
            x(bf>1)=1;
            x(bf>bfthresh)=2;
            x(bf<1/bfthresh)=-1;

            for z=-1:2
                plot(minmax(timevect),(st+lh*z+li*(i-1))*[1 1],'-','Color',zc,'LineWidth',ms);
            end
            plot(timevect,st+lh*x+li*(i-1),'o','Color',co(ii,:),'MarkerSize',ms);
            plot(timevect(x<0|x>1),st+lh*x(x<0|x>1)+li*(i-1),'o','Color',co(ii,:),'MarkerSize',ms,'MarkerFaceColor',co(ii,:));

            tt=text(max(timevect),st+li*(i-1),['  BF ' bn{i}],'VerticalAlignment','middle','FontSize',8,'Color',co(ii,:));
        end
        % BF diff
        for i=1:3
            ii=i;
            bf = bd{plotnr,i};
            x(bf>1)=1;
            x(bf>bfthresh)=2;
            x(bf<1/bfthresh)=-1;

            for z=-1:2
                plot(minmax(timevect),(st+lh*z+li*(i+2))*[1 1],'-','Color',zc,'LineWidth',ms);
            end
            plot(timevect,st+lh*x+li*(i+2),'o','Color',co2(ii,:),'MarkerSize',ms);
            plot(timevect(x<0|x>1),st+lh*x(x<0|x>1)+li*(i+2),'o','Color',co2(ii,:),'MarkerSize',ms,'MarkerFaceColor',co2(ii,:));

            tt=text(max(timevect),st+li*(i+2),['  BF ' bnd{i}],'VerticalAlignment','middle','FontSize',8,'Color',co2(ii,:));
        end
        legend(h,'Location','NE');
        
        a.XLim=(minmax(timevect));
        a.YLim=([.41 .6]);
        a.YTick=.5:.05:.6;
        a.XTick=0:100:1000;
        leg=legend();leg.Box='off';leg.Orientation='Vertical';
        xlabel('time (ms)')
        ylabel('accuracy')
        a.YLabel.Position = [-187 mean(a.YTick) -1];
        
        if plotnr==1
            % bf box
            a2 = axes('Position',[a.Position(1)-.095 a.Position(2) .07 .09],'box','on');
            %leg.Position = a2.Position+[0 a2.Position(4) 0 -.12];
            a2.XTick=[];a2.YTick=[];
            a2.YLim=[.05 1.1];a2.XLim=[0 1];hold on
            locs = linspace(0.1,.9,5);
            mlocs = movmean(locs,2);mlocs = mlocs(2:end);
            tt = fliplr({sprintf('BF>%i',bfthresh),'BF>1','BF<1',sprintf('BF<1/%i',bfthresh)});
            bfc = co(3,:);
            for y=1:length(locs)
                if y==3
                    plot([0.05 .95],locs(y)*[1 1],'-.','Color',bfc,'LineWidth',.8)
                else
                    plot([0.05 .95],locs(y)*[1 1],'Color',bfc,'LineWidth',.8)
                end
                if y<5
                    if ismember(y,[2 3])
                        plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6)
                    else
                        plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6,'MarkerFaceColor',bfc)
                    end
                    t=text(1/3,mlocs(y),tt{y},'Color',bfc,'FontSize',9);
                    
                end
            end
            t=text(.05,1,'Bayes Factor:','FontSize',10,'Color',bfc);
        end
        
    end
end

% bf box
% a2 = axes('Position',[.91 .84 .071 .1],'box','on');
% %leg.Position = a2.Position+[0 a2.Position(4) 0 -.12];
% a2.XTick=[];a2.YTick=[];
% a2.YLim=[.05 1.1];a2.XLim=[0 1];hold on
% locs = linspace(0.1,.9,5);
% mlocs = movmean(locs,2);mlocs = mlocs(2:end);
% tt = fliplr({sprintf('BF>%i',bfthresh),'BF>1','BF<1',sprintf('BF<1/%i',bfthresh)});
% bfc = 'k';
% for y=1:length(locs)
%     if y==3
%         plot([0.05 .95],locs(y)*[1 1],'-.','Color',bfc,'LineWidth',.8)
%     else
%         plot([0.05 .95],locs(y)*[1 1],'Color',bfc,'LineWidth',.8)
%     end
%     if y<5
%         if ismember(y,[2 3])
%             plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6)
%         else
%             plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6,'MarkerFaceColor',bfc)
%         end
%         t=text(1/3,mlocs(y),tt{y},'Color',bfc,'FontSize',9);
%     end
%     t=text(.1,1,'Bayes Factor:','FontSize',10);
% end

%%
fn = 'figures/figure_decoding_pairwise_half_sequence';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');


