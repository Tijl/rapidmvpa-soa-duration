%% plot_behaviour

%% load data
T = table();
for s=1:20
    load(sprintf('Archive/data_behavioural/data_subject_%i.mat',s),'data','p');
    
    st = [1 1+find(diff(data.sequencenumber))];
    S = struct();
    S.subject = s*ones(length(data.rt),1);
    S.sequencenumber = data.sequencenumber(st)';
    S.condition = data.condition(st)';
    S.querystim = data.querystim';
    S.response = data.response';
    S.correct = data.correct';
    S.rt = data.rt';
    T = [T; struct2table(S)];
    stimuli = p.stimuli;
end

%% get median RT and accuracy per person
prt=[];pchoice=[];

for i = 1:20
    for av= 1:2 % vehicle then animal
        idx = find(T.subject==i&ismember(T.querystim,(((av-1)*12+1):((av-1)*12+12))));
        
        prt(i,av) = median(T.rt(idx),'omitnan');
        pchoice(i,av) = mean(T.correct(idx),'omitnan');
    end
end
mean(pchoice)
std(pchoice)
mean(prt)

