%% plot image sequence timelines

%% load images
imf = dir('stimuli_resized/stim*');
ims={};alpha={};
rng(1);
order = randperm(24);
for i=1:24
    [ii,~,alpha{i}] = imread(fullfile(imf(order(i)).folder,imf(order(i)).name));
    if ndims(ii)<3
        ii = repmat(ii,1,1,3);
    end
    ims{i} = ii;
end

f=figure(1);clf
f.Position = [f.Position([1,2]) 1000 850];f.Resize='off';f.PaperPositionMode='auto';
drawnow

%bottom row: plot sequence
a=axes('Units','pixels','Position',[205 50 750 230]);hold on
a.YDir='reverse';
aw = .5;
axis equal
ff = 4;%factor decrease
a.XLim=[aw (96/ff)+aw];
a.YLim=[aw 5+aw];
a.XTick=[.5:1:(24/ff) (24/ff+.5):2:(48/ff) (48/ff+.5):4:(96/ff)];
a.TickDir='out';
a.XTickLabel=a.XTick*50-25;
a.XTickLabelRotation=45;
a.YTick=1:5;
a.FontSize=18;
a.TickLength=[.001 .001];
a.YTickLabel={'duration: 200 SOA: 200','duration: 100 SOA: 200',...
    'duration: 50 SOA: 200','duration: 50 SOA: 100',...
    'duration: 50 SOA: 50'};
a.XLabel.String='time (ms)';

onsets = [...
    1:4:96;...
    1:4:96;...
    1:4:96;...
    1:2:48;...
    1:24;...
];
offsets = [...
    4:4:96;...
    2:4:96;...
    1:4:96;...
    1:2:48;...
    1:24;...
];

for y=1:5
    m = 96/ff;
    if y==4;m=48/ff;end
    if y==5;m=24/ff;end
    fill([1-aw m+aw m+aw 1-aw],y+[-aw -aw +aw +aw],.75*[1 1 1],'LineWidth',1)
    for i=1:(24/ff)
        x = onsets(y,i);
        x2 = offsets(y,i);
        fill([x-aw x2+aw x2+aw x-aw],y+[-aw -aw aw aw],.9*[1 1 1],'LineWidth',1)
        image('XData',x+.9*[-aw aw],'Ydata',y+.9*[-aw aw],'CData',(ims{i}),'AlphaData',(alpha{i}))
        
    end
end

% Top row: plot all stimuli
imf = dir('stimuli_resized/stim*');
ims={};alpha={};
for i=1:24
    [ii,~,alpha{i}] = imread(fullfile(imf(i).folder,imf(i).name));
    if ndims(ii)<3
        ii = repmat(ii,1,1,3);
    end
    ims{i} = ii;
end

a=axes('Units','pixels','Position',[200 250 600 500]);hold on
a.YDir='reverse';
aw = .5;
axis equal
ff=4;
a.YLim=[aw 24/(24/ff)+aw];
a.XLim=[aw (24/ff)+aw];
axis off
for i=1:24
    y = 1+mod(i-1,ff);
    x = ceil(i/ff);
    image('XData',x+.9*[-aw aw],'Ydata',y+.9*[-aw aw],'CData',(ims{i}),'AlphaData',(alpha{i}))
end
co = viridis(4);
cnm={'birds','dogs','fish','boats','cars','planes'};
for i=1:6
    text(i,.3,cnm{i},'HorizontalAlignment','center','FontSize',20,'Color',co(2,:))
end
cnm={'animals','vehicles'};
for i=1:2
    text(2+3*(i-1),0,cnm{i},'HorizontalAlignment','center','FontSize',20,'Color',co(1,:))
end

cnm={'animacy','object','image'};
for i=1:3
    text(-.23,-.3+.3*i,cnm{i},'HorizontalAlignment','right','FontSize',20,'Color',co(i,:))
end

annotation('arrow','Units','pixels','Position',[140 747 60 0],'Color',co(1,:));
annotation('arrow','Units','pixels','Position',[140 718 60 0],'Color',co(2,:));
annotation('arrow','Units','pixels','Position',[140 685 60 -30],'Color',co(3,:));
annotation('arrow','Units','pixels','Position',[140 685 60 -80],'Color',co(3,:));
annotation('arrow','Units','pixels','Position',[140 685 60 -190],'Color',co(3,:));

t=annotation('textbox','Units','pixels','Position',[1 725 50 50],'String','A','FontSize',40,'LineStyle','none');
t=annotation('textbox','Units','pixels','Position',[1 250 50 50],'String','B','FontSize',40,'LineStyle','none');

%%
fn = 'figures/figure_design';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');

