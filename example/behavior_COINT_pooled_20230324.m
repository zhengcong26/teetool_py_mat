
clc;clear;close all;

currentFolder = pwd;
namelist = dir('**/*.mat');
J_tran_4_separate = [];
J_tran_4 = cell(1,4);
J_tran_group = cell(4,4);
J_velocity_list = [];
n=1;
% x=[];y=[];temp=[];
for ses = 1:length(namelist)
    if ~isempty(strfind(namelist(ses).name,'25msGaussian_4s')) %&& ~isempty(strfind(namelist(ses).name,'0815'))
        tic
        %         load(namelist(ses).name,'stim_list','J_GO_bin','J_MO_bin','J_TO_bin')
        load(namelist(ses).name,'stim_list','velocity_list','vicon')
        a = namelist(ses).name;
        sesname{n,1} = a(1:10);
        
        %% update velocity
        for i = 1:length(stim_list(:,1))
            % velocity 2
            [r1,~] = find(vicon(:,48) == stim_list(i,1) & vicon(:,50) == 8); % find peak index
            if ~isempty(r1)
                velocity_list{i,2}=[];
                velocity_list{i,2}(:,1) = vicon(r1-20:r1+20,51);
            end
            
            % velocity 4
            [r2,~] = find(vicon(:,48) == stim_list(i,1) & vicon(:,49) == 5); % find Go index
            if ~isempty(r2)
                velocity_list{i,4}=[];
                velocity_list{i,4}(:,1) = vicon(r2:r2+80,51);
            end
            
            [r3,~] = find(vicon(:,48) == stim_list(i,1) & vicon(:,50) == 61); % find MO index
            if ~isempty(r3)
                velocity_list{i,8}(:,1:3) = vicon(r3-10:r3+40,3:5);
            end
        end
        %% trialconstraint
        switch sesname{n,1}
            case 'G_20220905'
                trialconstraint = 582;
            case 'G_20220906'
                trialconstraint = 779;
            case 'G_20220907'
                trialconstraint = 642;
            case 'G_20220908'
                trialconstraint = 648;
            case 'G_20220910'
                trialconstraint = 646;
            case 'G_20220912'
                trialconstraint = 613;
            case 'G_20220914'
                trialconstraint = 675;
            case 'G_20220915'
                trialconstraint = 604;
            case 'G_20220916'
                trialconstraint = 610;
            case 'G_20220917'
                trialconstraint = 662;
            case 'G_20221005'
                trialconstraint = 614;
            case 'G_20221007'
                trialconstraint = 580;
            case 'G_20221009'
                trialconstraint = 792;
            case 'G_20221012'
                trialconstraint = 597;
            case 'L_20220809'
                trialconstraint = 848;
            case 'L_20220810'
                trialconstraint = 634;
            case 'L_20220812'
                trialconstraint = 1000;
            case 'L_20220813'
                trialconstraint = 666;
            case 'L_20220815'
                trialconstraint = 675;
            case 'L_20220818'
                trialconstraint = 644;
            case 'L_20220819'
                trialconstraint = 674;
            case 'L_20220821'
                trialconstraint = 848;
            case 'L_20220824'
                trialconstraint = 611;
            case 'L_20220827'
                trialconstraint = 747;
            case 'L_20220828'
                trialconstraint = 657;
            case 'L_20220909'
                trialconstraint = 1000;
            case 'L_20220912'
                trialconstraint = 1200;
            case 'L_20220913'
                trialconstraint = 1000;
            case 'L_20220914'
                trialconstraint = 1000;
            case 'L_20220915'
                trialconstraint = 895;
            case 'L_20220916'
                trialconstraint = 756;
            case 'G_20220930'
                trialconstraint = 1000; % no ICMS use all trials
            case 'G_20221001'
                trialconstraint = 1000; % no ICMS use all trials
            case 'L_20220822'
                trialconstraint = 1000; % no ICMS use all trials
            case 'L_20220823'
                trialconstraint = 1000; % no ICMS use all trials
        end
        %% CO INT paired
        stim_list(stim_list(:,4)>230&stim_list(:,4)<240,4) = 216;
        stim_list(stim_list(:,4)>930&stim_list(:,4)<940,4) = 916;
        
        stim_list_select = sortrows(stim_list(stim_list(:,4)>0&stim_list(:,5)==0&stim_list(:,1)<=trialconstraint,:),1);
        
        CO = stim_list_select(stim_list_select(:,3)==0,:);
        INT = cell(size(CO,1),1);
        
        for s = 1:size(CO,1)
            condition_non = stim_list_select(stim_list_select(:,17)>CO(s,17)-20 & stim_list_select(:,17)<CO(s,17)+20 &...
                stim_list_select(:,3)~=0 & stim_list_select(:,4)>CO(s,4)-100 & stim_list_select(:,4)<CO(s,4)+100,:);
            INT{s,1} = [INT{s,1};condition_non];
        end
        
        INT_mat = cell2mat(INT);

        [~,ia,~] = unique(INT_mat(:,1),'rows');
        INT_unique = INT_mat(ia,:);
        INT_unique = sortrows(INT_unique,1);
        %% tran_group
        tran = [CO;INT_unique];
        
        hand_loc=zeros(size(tran,1),1);
        hand_loc(0<tran(:,17) & tran(:,17)<90)=1;
        hand_loc(90<tran(:,17) & tran(:,17)<180)=2;
        hand_loc(180<tran(:,17) & tran(:,17)<270)=3;
        hand_loc(270<tran(:,17) & tran(:,17)<360)=4;
        hand_loc(360<tran(:,17))=1;
        
        sp=unique(tran(:,3));
        movingspeed = setdiff(sp,0);
        sp=[0,movingspeed];
        tran_group_1 = [];tran_group_2 = [];tran_group_3 = [];
        for i = 1:2
            for j = 1:4
                tran_group_1{j,i}=tran(tran(:,4)<500&tran(:,3)==sp(i)&hand_loc==j,:);% short
                tran_group_2{j,i}=tran(tran(:,4)>500&tran(:,3)==sp(i)&hand_loc==j,:);% long
                tran_group_3{j,i}=tran(tran(:,3)==sp(i)&hand_loc==j,:);
            end
        end
        
        tran_group = [tran_group_1,tran_group_2]; % row:4Q, column:short CO,short INT,long CO,long INT
        COINT = tran_group_3;
        %% for graphpad prism
        
        tran_4=cell(1,4); % short CO,short INT,long CO,long INT
        for i = 1:4
            if ~isempty(tran_group(:,i))
                tran_4{1,i} = cell2mat(tran_group(:,i));
            else
                tran_4{1,i} = nan(1,72);
            end
        end
        
        J_tran_4_separate = [J_tran_4_separate;tran_4];
        J_tran_4 = cellfun(@(x,y) ([x;y]),J_tran_4,tran_4,'UniformOutput',0);
        J_tran_group = cellfun(@(x,y) ([x;y]),J_tran_group,tran_group,'UniformOutput',0);
        J_velocity_list = [J_velocity_list;velocity_list];
        n=n+1;
        toc
    end
end

%% para
para = 13; % RT
% para = 48; % MT
% para = 49; % RT+MT
% para = 36; % V
% para = 30; % angle
tran_para = cellfun(@(x) (x(:,para)),J_tran_4,'UniformOutput',0);
tran_group_para = cellfun(@(x) (x(:,para)),J_tran_group,'UniformOutput',0);

% para = 18;
% tran_para = cellfun(@(x) (abs(x(:,para)*10.2/16)),J_tran_4,'UniformOutput',0);
% tran_group_para = cellfun(@(x) (abs(x(:,para)*10.2/16)),J_tran_group,'UniformOutput',0);

%% statistic
tran_para4 = nan(2000,4);
for i=1:4
    tran_para4(1:size(tran_para{1,i},1),i) = tran_para{1,i};
end

n = 1;
tran_para4_group = nan(500,16);
for i=1:4
    for j = 1:4
        tran_para4_group(1:size(tran_group_para{i,j},1),n) = tran_group_para{i,j};
        n = n+1;
    end
end

x=tran_para4;
behavelist_4=[];diff_p_4=[];
for i = 1:2:3
    diff_p_4(1,i)=median(x(:,i+1),'omitnan')-median(x(:,i),'omitnan');
    diff_p_4(1,i+1)=ranksum(x(:,i),x(:,i+1));
end
for j = 1:4
    behavelist_4(1,j) = length(x(~isnan(x(:,j)),j));
    behavelist_4(3,j) = median(x(:,j),'omitnan');
    temp = bootci(1000,@median,x(~isnan(x(:,j)),j));
    behavelist_4(2,j) = temp(1);
    behavelist_4(4,j) = temp(2);
end
behavelist_4(5,:) = diff_p_4(1,:);

x = tran_para4_group;
diff_p=[];behavelist=[];
for i = 1:2:15
    diff_p(1,i)=median(x(:,i+1),'omitnan')-median(x(:,i),'omitnan');
    diff_p(1,i+1)=ranksum(x(:,i),x(:,i+1));
end
for j = 1:16
    behavelist(1,j) = length(x(~isnan(x(:,j)),j));
    behavelist(3,j) = median(x(:,j),'omitnan');
    temp = bootci(1000,@median,x(~isnan(x(:,j)),j));
    behavelist(2,j) = temp(1);
    behavelist(4,j) = temp(2);
end
behavelist(5,:) = diff_p(1,:);

fortable=[behavelist_4,behavelist];

% %% statistic
% x=tran_para{1,1}; % CO
% y=tran_para{1,2}; % INT
% [p,h,stats]=ranksum(x,y)
% [h,p] = kstest2(x,y)
% % [p,h,stats]=ranksum(x,y,'method','exact');
% [p,h,dif] = bootstrap2(x,y)
% 
% 
% x=tran_para{1,3}; % CO
% y=tran_para{1,4}; % INT
% [p,h,~]=ranksum(x,y)
% [h,p] = kstest2(x,y)
% % [p,h,stats]=ranksum(x,y,'method','exact');
% [p,h,dif] = bootstrap2(x,y)


%% reciprobit regress 20230417

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linest={'--','-'};

y1 = tran_para{1,1};y2 = tran_para{1,2};y3 = tran_para{1,3};y4 = tran_para{1,4};
betas=[];p=[];r2=[];
figure
set(gcf,'unit','centimeters','position',[10 5 20 8]);
hold on
% 1
[betas(:,1), xa,ya] = reciprobit_noplot(y1,[],[],'range',[50 inf]);
[h,p{1,1},r2(:,1)] = reciprobit_regress(xa,ya);
h(1).MarkerEdgeColor = colo(1,:);
h(1).MarkerSize = 1.5;
% h(2).Color = colo(1,:);
% h(2).LineStyle = '--';
% h(2).LineWidth = 1;

h(3).Color = colo(1,:);
h(3).LineStyle = linest{1,1};
h(3).LineWidth = 1;
h(2).Color = colo(1,:);
h(2).LineStyle = linest{1,1};
h(2).LineWidth = 0.1;

% 2
[betas(:,2), xa,ya] = reciprobit_noplot(y2,[],[],'range',[50 inf]);
[h,p{1,2},r2(:,2)] = reciprobit_regress(xa,ya);
h(1).MarkerEdgeColor = colo(2,:);
h(1).MarkerSize = 1.5;
h(2).Color = colo(2,:);
h(2).LineStyle = linest{1,1};
h(2).LineWidth = 1;

% 3
[betas(:,3),xa,ya] = reciprobit_noplot(y3,[],[],'range',[50 inf]);
[h,p{1,3},r2(:,3)] = reciprobit_regress(xa,ya);
h(1).MarkerEdgeColor = colo(1,:);
h(1).MarkerSize = 1.5;
h(3).Color = colo(1,:);
h(3).LineStyle = '-';
h(3).LineWidth = 1;
h(2).Color = colo(1,:);
h(2).LineStyle = '-';
h(2).LineWidth = 0.1;

% 4
[betas(:,4), xa,ya] = reciprobit_noplot(y4,[],[],'range',[50 inf]);
[h,p{1,4},r2(:,4)] = reciprobit_regress(xa,ya);
h(1).MarkerEdgeColor = colo(2,:);
h(1).MarkerSize = 1.5;
% h(3).Color = colo(2,:);
% h(3).LineStyle = '-';
% h(3).LineWidth = 1;
% h(2).Color = colo(2,:);
% h(2).LineStyle = '-';
% h(2).LineWidth = 0.1;
h(2).Color = colo(2,:);
h(2).LineStyle = '-';
h(2).LineWidth = 1;

xlims=xlim(); 
xlim([xlims(1),0]);           % ensure the RT axis goes up to infinity
xvals=get(gca,'XTick'); 
xvals(xvals==0)=-0; % replace 0 by -0
% set up correct RT labels on the x-axis (inverse scale -1/RT)
set(gca,'XTickLabel', cellfun( @(x)sprintf('%.3g',-1./x), num2cell(xvals) , 'UniformOutput', 0));
xlabel('RT'); 
myerf=@(x)(erf(x)+1).*50; 
myerfinv=@(x)erfinv(x./50-1);
yvals=linspace(myerf(-3),myerf(3),length(get(gca,'YTick')));
set(gca,'YTick', myerfinv(yvals) );         % and use probit proportions
set(gca,'YTickLabel',cellfun( @(x)sprintf('%.3g',x), num2cell(yvals), 'UniformOutput',0));
ylabel('cumulative probability (%)');       % add a dotted line at P = 50%
h=line([xlims(1),0],[0,0]); set(h,'LineStyle',':','LineWidth',1, 'Color',get(gca,'XColor'));
ylim([-4,4])

% 
filename = 'L -120';
print(1, '-dpng', '-r600', [datestr(now,29) '_reciprobit_regress_' filename])
close all
%% reciprobit regress statistic 20230417
p=[];r2=[];
for rep = 1:100
    y1 = datasample(tran_para{1,1},size(tran_para{1,1},1));
    y2 = datasample(tran_para{1,2},size(tran_para{1,2},1));
    y3 = datasample(tran_para{1,3},size(tran_para{1,3},1));
    y4 = datasample(tran_para{1,4},size(tran_para{1,4},1));
    
    [~, xa,ya] = reciprobit_noplot(y1,[],[],'range',[50 inf]);
    [~,p{rep,1},r2(:,1)] = reciprobit_regress(xa,ya);
    close all
    % 2
    [~, xa,ya] = reciprobit_noplot(y2,[],[],'range',[50 inf]);
    [~,p{rep,2},r2(:,2)] = reciprobit_regress(xa,ya);
    close all
    % 3
    [~,xa,ya] = reciprobit_noplot(y3,[],[],'range',[50 inf]);
    [~,p{rep,3},r2(:,3)] = reciprobit_regress(xa,ya);
    close all
    % 4
    [~, xa,ya] = reciprobit_noplot(y4,[],[],'range',[50 inf]);
    [~,p{rep,4},r2(:,4)] = reciprobit_regress(xa,ya);
    close all
end

slope = cell2mat(cellfun(@(x) (x(2,1)),p,'UniformOutput',0));
intercept = cell2mat(cellfun(@(x) (x(2,2)),p,'UniformOutput',0));
%% reciprobit figure 20230417 
% colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
% linest={'--','-'};
% y1 = tran_para{1,1};y2 = tran_para{1,2};y3 = tran_para{1,3};y4 = tran_para{1,4};
% 
% figure
% set(gcf,'unit','centimeters','position',[10 5 7 5]);
% hold on
% % 1
% [betas(:,1), graphhandles, ~,~] = reciprobit(y1,[],[],'range',[50 inf]);
% for i = 1:2
%     graphhandles(i).Color = colo(1,:);
% end
% graphhandles(1).Marker = 'o';graphhandles(1).MarkerSize = 0.1;graphhandles(2).LineStyle = '--';graphhandles(2).LineWidth = 1;
% 
% % 2
% [betas(:,2), graphhandles, ~,~] = reciprobit(y2,[],[],'range',[50 inf]);
% for i = 1:2
%     graphhandles(i).Color = colo(2,:);
% end
% graphhandles(1).Marker = 'o';graphhandles(1).MarkerSize = 0.2;graphhandles(2).LineStyle = '--';graphhandles(2).LineWidth = 1;
% 
% % 3
% [betas(:,3), graphhandles,~,~] = reciprobit(y3,[],[],'range',[50 inf]);
% for i = 1:2
%     graphhandles(i).Color = colo(1,:);
% end
% graphhandles(1).Marker = '+';graphhandles(1).MarkerSize = 0.3;graphhandles(2).LineStyle = '-';graphhandles(2).LineWidth = 1;
% 
% [betas(:,4), graphhandles, x_all,y_all] = reciprobit(y4,[],[],'range',[50 inf]);
% for i = 1:2
%     graphhandles(i).Color = colo(2,:);
% end
% graphhandles(1).Marker = '+';graphhandles(1).MarkerSize = 0.3;graphhandles(2).LineStyle = '-';graphhandles(2).LineWidth = 1;
% 
% xlims=xlim(); 
% xlim([xlims(1),0]);           % ensure the RT axis goes up to infinity
% xvals=get(gca,'XTick'); 
% xvals(xvals==0)=-0; % replace 0 by -0
% set(gca,'XTickLabel', cellfun( @(x)sprintf('%.3g',-1./x), num2cell(xvals) , 'UniformOutput', 0));
% ylim([-4 4])
% 
% filename = 'L -120';
% print(1, '-dpng', '-r600', [datestr(now,29) '_reciprobit_' filename])
% close all
%% reciprobit statistic 20230417 
% betas=[];
% for rep = 1:100
%     y1 = datasample(tran_para{1,1},size(tran_para{1,1},1));
%     y2 = datasample(tran_para{1,2},size(tran_para{1,2},1));
%     y3 = datasample(tran_para{1,3},size(tran_para{1,3},1));
%     y4 = datasample(tran_para{1,4},size(tran_para{1,4},1));
%     
%     [betas{1,1}(:,rep), ~, ~,~] = reciprobit(y1,[],[],'range',[50 inf]);
%     close all
%     
%     [betas{1,2}(:,rep), ~, ~,~] = reciprobit(y2,[],[],'range',[50 inf]);
%     close all
%     
%     [betas{1,3}(:,rep), ~,~,~] = reciprobit(y3,[],[],'range',[50 inf]);
%     close all
%     
%     [betas{1,4}(:,rep), ~, ~,~] = reciprobit(y4,[],[],'range',[50 inf]);
%     close all
%     
% end
% 
% intercept = cell2mat(cellfun(@(x) (transpose(x(1,:))),betas,'UniformOutput',0));
% slope = cell2mat(cellfun(@(x) (transpose(x(2,:))),betas,'UniformOutput',0));


%% cdf all Q

filename = 'L -120';
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

figure
ax = axes;
x = tran_para{1,1};
h1=cdfplot_my(x);
set(h1,'LineStyle', '--', 'Color', colo(1,:),'LineWidth',1) % CO short
hold on
scatter(median(x),0.05,20,'o','MarkerEdgeColor',colo(1,:))
x = tran_para{1,2};
h2=cdfplot_my(x);
set(h2,'LineStyle', '--', 'Color',colo(2,:),'LineWidth',1) % INT short
scatter(median(x),0.05,20,'o','MarkerEdgeColor',colo(2,:))
x = tran_para{1,3};
h3=cdfplot_my(x);
set(h3,'LineStyle', '-', 'Color', colo(3,:),'LineWidth',1) % CO long
scatter(median(x),0.05,20,'o','filled','MarkerEdgeColor',colo(3,:),'MarkerFaceColor',colo(3,:))
x = tran_para{1,4};
h4=cdfplot_my(x);
set(h4,'LineStyle', '-', 'Color',colo(4,:),'LineWidth',1) % INT long
scatter(median(x),0.05,20,'o','filled','MarkerEdgeColor',colo(4,:),'MarkerFaceColor',colo(4,:))
xlim([0 600])
% xlim([200 250])
box off
grid off
ylabel('cdf','FontName','Arial');
xlabel('reaction time (ms)','FontName','Arial');
ax.TickDir='out';
set(gcf,'unit','centimeters','position',[10 5 7 5]);
set(gca,'Position',[.2 .2 .7 .7]);
offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'

print(1, '-dpng', '-r600', [datestr(now,29) '_COINT_RT_cdf_' filename])
close all
%% cdf 4Q
for i = 1:4
    figure
    ax = axes;
    x = tran_group_para{i,1};
    h1=cdfplot_my(x);
    set(h1,'LineStyle', '--', 'Color', [0.06275	0.30588	0.5451],'LineWidth',1)
    hold on
    x = tran_group_para{i,2};
    h2=cdfplot_my(x);
    set(h2,'LineStyle', '--', 'Color',[0.80392	0.33333	0.33333],'LineWidth',1)
    x = tran_group_para{i,3};
    h3=cdfplot_my(x);
    set(h3,'LineStyle', '-', 'Color', [0.06275	0.30588	0.5451],'LineWidth',1)
    hold on
    x = tran_group_para{i,4};
    h4=cdfplot_my(x);
    set(h4,'LineStyle', '-', 'Color',[0.80392	0.33333	0.33333],'LineWidth',1)
    xlim([0 600])
    box off
    grid off
    ylabel('cdf','FontName','Arial');
    xlabel('reaction time (ms)','FontName','Arial');
    ax.TickDir='out';
    set(gcf,'unit','centimeters','position',[10 5 7 5]);
    set(gca,'Position',[.2 .2 .7 .7]);
    offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    
    print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_' num2str(i)])
    close all
end

%% cdf larger
filename = 'L -120';
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

figure
ax = axes;
x = tran_para{1,1};
h1=cdfplot_my(x);
set(h1,'LineStyle', '--', 'Color', colo(1,:),'LineWidth',1) % CO short
hold on
scatter(median(x),0.05,20,'o','MarkerEdgeColor',colo(1,:))
x = tran_para{1,2};
h2=cdfplot_my(x);
set(h2,'LineStyle', '--', 'Color',colo(2,:),'LineWidth',1) % INT short
scatter(median(x),0.05,20,'o','MarkerEdgeColor',colo(2,:))
x = tran_para{1,3};
h3=cdfplot_my(x);
set(h3,'LineStyle', '-', 'Color', colo(3,:),'LineWidth',1) % CO long
scatter(median(x),0.05,20,'o','filled','MarkerEdgeColor',colo(3,:),'MarkerFaceColor',colo(3,:))
x = tran_para{1,4};
h4=cdfplot_my(x);
set(h4,'LineStyle', '-', 'Color',colo(4,:),'LineWidth',1) % INT long
scatter(median(x),0.05,20,'o','filled','MarkerEdgeColor',colo(4,:),'MarkerFaceColor',colo(4,:))
xlim([200 250])
% xlim([240 280])
ylim([0 1])
box on
grid off
% ylabel('cdf','FontName','Arial');
% xlabel('reaction time (ms)','FontName','Arial');
ax.TickDir='out';
set(gcf,'unit','centimeters','position',[10 5 2 3]);
set(gca,'Position',[.2 .2 .7 .7]);
% offsetaxis(ax, 'y', 0.01); % Default, y-offset = 0.1 and yloc = 'l'
% offsetaxis(ax, 'x', 0.01); % Default, y-offset = 0.1 and yloc = 'l'

print(1, '-dpng', '-r600', [datestr(now,29) '_COINT_RT_cdf_larger_' filename])
close all

%% histogram RT 20230717
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
style = {'--','--','-','-'};

figure
for i = 1:4
    x = tran_para4(:,i);
    histogram_time(25,x,colo(i,:),style{i},2)
    hold on
end




%% figure single trial trajectory no initial position deletion
figure
set(gcf,'position',[30,300,410,400]);
hold on
axis square
axis off
plot([-50 0],[60 60],'k','LineWidth',1)
text(-40,52,'5 cm','FontSize',14)
%% COINT
for q = 1:4
    
    %% CO
    n=1;xC=[];yC=[];
    for k = 1:size(COINT{q,1},1)
        for i = 1:length(velocity_list)
            if ~isempty(velocity_list{i,1})
                if velocity_list{i,1}(1,3) == COINT{q,1}(k,1) %&& velocity_list{i,1}(1,4) == COINT{q,1}(k,32)
                    try
                        xC(:,n) = -velocity_list{i,3}(1:10,1);
                        yC(:,n) = velocity_list{i,3}(1:10,2);
                        n = n+1;
                    catch
                    end
                    hold on
                end
            end
        end
    end
    
    xC_m = mean(xC,2,'omitnan'); % mean
    yC_m = mean(yC,2,'omitnan'); % mean
    
    %% plot
    plot(xC_m,yC_m,'color',[0.06275	0.30588	0.5451],'linewidth',3);
    
    for i = 1:size(xC,2)
        p = plot(smoothdata(xC(:,i),'gaussian',5), smoothdata(yC(:,i),'gaussian',5),'color',[0.06275	0.30588	0.5451],'linewidth',1.2);
        p.Color(4) = 0.25;
        hold on
    end
    
    %% ST
    n=1;xS=[];yS=[];
    for k = 1:size(COINT{q,2},1)
        for i = 1:length(velocity_list)
            if ~isempty(velocity_list{i,1})
                if velocity_list{i,1}(1,3) == COINT{q,2}(k,1) %&& velocity_list{i,1}(1,4) == COINT{q,2}(k,32)
                    try
                        xS(:,n) = -velocity_list{i,3}(1:10,1);
                        yS(:,n) = velocity_list{i,3}(1:10,2);
                        n = n+1;
                    catch
                    end
                    hold on
                end
            end
        end
    end
    
    xS_m = mean(xS,2,'omitnan'); % mean
    yS_m = mean(yS,2,'omitnan'); % mean
    
    %% plot
    plot(xS_m,yS_m,'color',[0.80392	0.33333	0.33333],'linewidth',3);
    
    for i = 1:size(xS,2)
        p = plot(smoothdata(xS(:,i),'gaussian',5), smoothdata(yS(:,i),'gaussian',5),'color',[0.80392	0.33333	0.33333],'linewidth',1.2);
        p.Color(4) = 0.25;
        hold on
    end
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_trj_no_rotation'])
close all

%% figure mean CI trajectory
%% COINT
XYr = [];XYr_trial = [];
for q = 1:4
    for c = 1:4
        
        n = 1;tjC = [];
        for k = 1:size(J_tran_group{q,c},1)
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i,1})
                    if J_velocity_list{i,1}(1,3) == J_tran_group{q,c}(k,1) && J_velocity_list{i,1}(1,4) == J_tran_group{q,c}(k,32)
                        try
                            x0 = -J_velocity_list{i,3}(1,1);
                            y0 = J_velocity_list{i,3}(1,2);
                            tjC(1:10,n) = -J_velocity_list{i,3}(1:10,1)-x0;
                            tjC(11:20,n) = J_velocity_list{i,3}(1:10,2)-y0;
                            n = n+1;
                            %                                         plot(-velocity_list{i,3}(1:10,1)-x0, velocity_list{i,3}(1:10,2)-y0,'k','linewidth',1.5) ;
                        catch
                        end
                    end
                end
            end
        end
        
        xC = mean(tjC(1:10,:),2,'omitnan'); % mean
        yC = mean(tjC(11:end,:),2,'omitnan');
        %% rotate
        [XYr{q,c},~] = rotatetrj(xC,yC);
        XYr{q,c} = [smoothdata(XYr{q,c}(1,:)','gaussian',4),smoothdata(XYr{q,c}(2,:)','gaussian',4)];
        
        for i = 1:size(tjC,2)
            [XYr_trial_temp,~] = rotatetrj(tjC(1:10,i),tjC(11:20,i));
            XYr_trial{q,c}{1,1}(:,i) = smoothdata(XYr_trial_temp(1,:)','gaussian',4);
            XYr_trial{q,c}{1,2}(:,i) = smoothdata(XYr_trial_temp(2,:)','gaussian',4);
        end
               
    end
end
%% plot rotate
figure
set(gcf,'position',[30,300,620,230]);
pos = [2,1,3,4];
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for i = 1:4
    subplot(2,2,pos(i))
    hold on
    for j = 1:4
        
        %         if i == 1
        %             M=[cos(0),-sin(0);sin(0),cos(0)];
        %         elseif i == 2 || i == 3 || i == 4
        %             M=[cos(2*pi),sin(2*pi);-sin(2*pi),cos(2*pi)];
        %         end
        %         R1=[XYr{i,j}(:,1),XYr{i,j}(:,2)];
        %         R2=M*R1';%旋转后坐标
        %         plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
        plot(XYr{i,j}(:,1),XYr{i,j}(:,2),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
%                 % single
%                 for k = 1:size(XYr_trial{i,j}{1,1},2)
%                     p = plot(XYr_trial{i,j}{1,1}(:,k),XYr_trial{i,j}{1,2}(:,k),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',1);
%                     if j == 1 || j == 3
%                         p.Color(4) = 0.8;
%                     else
%                         p.Color(4) = 0.4;
%                     end
%                 end
        
        %
        
    end
    y0=ylim;
    plot([-10 -10],[y0(1)+1 y0(1)+6],'k','LineWidth',1)

    if i == 1
        plot([-10 10],[y0(1)+1 y0(1)+1],'k','LineWidth',1)
    end
    %     axis equal
    axis off
    axis([-10 120 -inf inf])
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_trj_rotation'])
close all

%% trj plot rotate to origin position
figure
set(gcf,'position',[30,300,500,400]);
hold on
axis([-120 120 -120 120])
axis square
axis equal
axis off
plot([-50 0],[-70 -70],'k','LineWidth',1)
text(-40,-80,'5 cm','FontSize',14)
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for i = 1:4
    for j = 1:4
        if i == 1
            M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
        elseif i == 2
            M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
        elseif i == 3
            M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
        elseif i == 4
            M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
        end
        
                R1=[XYr{i,j}(:,1),XYr{i,j}(:,2)];
                R2=M*R1';%旋转后坐标
                plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
%         % single trials
%         for k = 1:size(XYr_trial{i,j}{1,1},2)
%             R1=[XYr_trial{i,j}{1,1}(:,k),XYr_trial{i,j}{1,2}(:,k)];
%             R2=M*R1';%旋转后坐标
%             p = plot(R2(1,:),R2(2,:),'Color',colo(j,:),'LineStyle',linestyle{j},'LineWidth',1);
%             if j == 1 || j == 3
%                 p.Color(4) = 0.8;
%             else
%                 p.Color(4) = 0.4;
%             end
%         end
        
    end
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_trj'])
close all

%% trj plot rotate to origin position single trial, also output rotrj for teetool 20230826
figure
set(gcf,'position',[30,300,500,400]);
hold on
axis([-120 120 -120 120])
axis square
axis equal
axis off
plot([-50 0],[-70 -70],'k','LineWidth',1)
text(-40,-80,'5 cm','FontSize',14)
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
rotrj = [];
for i = 1:4
    for j = 1:4
        if i == 1
            M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
        elseif i == 2
            M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
        elseif i == 3
            M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
        elseif i == 4
            M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
        end
        
%                 R1=[XYr{i,j}(:,1),XYr{i,j}(:,2)];
%                 R2=M*R1';%旋转后坐标
%                 plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
        % single trials
        for k = 1:size(XYr_trial{i,j}{1,1},2)
            R1=[XYr_trial{i,j}{1,1}(:,k),XYr_trial{i,j}{1,2}(:,k)];
            R2=M*R1';%旋转后坐标
            p = plot(R2(1,:),R2(2,:),'Color',colo(j,:),'LineStyle',linestyle{j},'LineWidth',1);
            rotrj{i,j}{k,1}=R2;
            if j == 1 || j == 3
                p.Color(4) = 0.8;
            else
                p.Color(4) = 0.4;
            end
        end
        
    end
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_single_trj'])
close all

%% try CI
XYr_trial_R=[];
XYr_trial_CI=[];
for i = 1:4
    for j = 1:4
        if i == 1
            M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
        elseif i == 2
            M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
        elseif i == 3
            M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
        elseif i == 4
            M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
        end
        
        for k = 1:size(XYr_trial{i,j}{1,1},2)
            R1=[XYr_trial{i,j}{1,1}(:,k),XYr_trial{i,j}{1,2}(:,k)];
            R2=M*R1';%旋转后坐标
            XYr_trial_R{i,j}{1,1}(:,k)=R2(1,:)';
            XYr_trial_R{i,j}{1,2}(:,k)=R2(2,:)';
        end
        
        x = XYr_trial_R{i,j}{1,1};
        y = XYr_trial_R{i,j}{1,2};
        x = x(:,all(~isnan(x)));
        y = y(:,all(~isnan(y)));
        
        for k = 1:10
            XYr_trial_CI{i,j}{1,1}(:,k) = bootci(1000,@mean,x(k,:));
            XYr_trial_CI{i,j}{1,2}(:,k) = bootci(1000,@mean,y(k,:));
        end
        
    end
end
%% try CI plot
figure
set(gcf,'position',[30,300,500,400]);
hold on
axis([-120 120 -120 120])
axis square
axis equal
axis off
plot([-50 0],[-70 -70],'k','LineWidth',1)
text(-40,-80,'5 cm','FontSize',14)
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for i = 1:4
    for j = 1:4
        if i == 1
            M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
        elseif i == 2
            M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
        elseif i == 3
            M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
        elseif i == 4
            M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
        end
        
        R1=[XYr{i,j}(:,1),XYr{i,j}(:,2)];
        R2=M*R1';%旋转后坐标
        plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
        for k = 1:size(XYr_trial_CI{i,j}{1,1},1)
            x = XYr_trial_CI{i,j}{1,1}(k,:);
            y = XYr_trial_CI{i,j}{1,2}(k,:);
            p = plot(x,y,'Color',colo(j,:),'LineStyle',linestyle{j},'LineWidth',1);
            p.Color(4) = 0.6;
        end
        
    end
end

% plot
%     plot(xC, yC,'b','linewidth',2) ;
%     hold on
%
%     CI = [];
%     for i = 1:size(tjC,1)
%         CI(i,:) = bootci(1000,@mean,tjC(i,~isnan(tjC(i,:))));
%     end
%
%     CI_up = [CI(1:10,1),CI(11:20,1)];
%     CI_down = [CI(1:10,2),CI(11:20,2)];
%
%     h=fill([CI_up(:,1);CI_down(end:-1:1,1)],[CI_up(:,2);CI_down(end:-1:1,2)],'b');
%     set(h,'FaceColor','b','FaceAlpha',0.3,'EdgeColor','none');




%% figure v GO mean
tj = [];
tran_v = J_tran_group;

figure
set(gcf,'position',[30,300,620,620]);
pos = [2,1,3,4];
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for a = 1:4
    subplot(2,2,pos(a))
    hold on
    for b = 1:4
        n = 1;
        tj = [];
        for k = 1:size(tran_v{a,b})
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i,1})
                    try
                        if J_velocity_list{i,1}(1,3) == tran_v{a,b}(k,1) && J_velocity_list{i,1}(1,4) == tran_v{a,b}(k,32)
                            tj(:,n) = J_velocity_list{i,4}(:,1);
                            n = n+1;
                        end
                    catch
                    end
                end
            end
        end
        
        y0 = mean(tj,2,'omitnan');
        y1 = y0';
        smoo = smoothdata(y1,'gaussian',8);
        t = 0:10:(length(smoo)-1)*10;
        
        plot(t,smoo,'color',colo(b,:),'LineStyle',linestyle{b},'LineWidth',1); % black
        
        std0 = std(tj,0,2,'omitnan')/sqrt(length(tj(1,:))-1); % std
        std1 = [y0-std0, y0(end:-1:1)+std0(end:-1:1)];
        std2(1,1:length(smoo)) = std1(1:length(smoo),1)';
        std2(1,length(smoo)+1:2*length(smoo)) = std1(1:length(smoo),2)';
        std2 = smoothdata(std2,'gaussian',8);
        
        if std0 ~= 0
            h=fill([t, t(end:-1:1)],std2,'b');
            set(h,'FaceColor',colo(b,:),'FaceAlpha',0.2,'EdgeColor','none');
        end
        
        axis([0 500 0 100])
        axis off
        if a == 1
            plot([0 0],[20 30],'k','LineWidth',1)
            plot([0 50],[20 20],'k','LineWidth',1)
            
            xlabel('Time from GO (ms)')
            ylabel('Hand speed (cm/s)')
        end
        
    end
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_velocity_GO'])
close all


%% figure v peak mean
tran_v = J_tran_group;
t = -200:10:200;

figure
% set(gcf,'position',[30,100,620,620]);
% set(gcf,'unit','centimeters','position',[10 5 10 10]);
set(gcf,'unit','centimeters','position',[10 5 12 10]);
pos = [2,1,3,4];
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for a = 1:4
    subplot(2,2,pos(a))
    hold on
    for b = 1:4
        n = 1;
        tj = [];
        for k = 1:size(tran_v{a,b})
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i,1})
                    try
                        if J_velocity_list{i,1}(1,3) == tran_v{a,b}(k,1) && J_velocity_list{i,1}(1,4) == tran_v{a,b}(k,32)
                            tj(:,n) = J_velocity_list{i,2}(:,1);
                            n = n+1;
                        end
                    catch
                    end
                end
            end
        end
        
        y0 = mean(tj,2,'omitnan');
        y1 = y0';
        smoo = smoothdata(y1,'gaussian',8);
        %         t = 0:10:(length(smoo)-1)*10;
        
        plot(t,smoo,'color',colo(b,:),'LineStyle',linestyle{b},'LineWidth',1.5); % black
        
        std2 = [];
        std0 = std(tj,0,2,'omitnan')/sqrt(length(tj(1,:))-1); % std
        std1 = [y0-std0, y0(end:-1:1)+std0(end:-1:1)];
        std2(1,1:length(smoo)) = std1(1:length(smoo),1)';
        std2(1,length(smoo)+1:2*length(smoo)) = std1(1:length(smoo),2)';
        std2 = smoothdata(std2,'gaussian',8);
        
        if std0 ~= 0
            h=fill([t, t(end:-1:1)],std2,'b');
            set(h,'FaceColor',colo(b,:),'FaceAlpha',0.2,'EdgeColor','none');
        end
        
        axis([-80 80 40 100])
%         axis off
        if a == 3
            plot([-200 -200],[50 75],'k','LineWidth',1)
            plot([-200 -140],[50 50],'k','LineWidth',1)
            
            xlabel('Time to Peak velocity (ms)')
            ylabel('Hand speed (cm/s)')
        end
        
    end
end

% print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_velocity_Peak'])
% close all


%% figure touchpoint rotate
figure
set(gcf,'unit','centimeters','position',[10 5 10 10]);
% set(gcf,'position',[30,300,500,500]);
pos = [2,1,3,4];
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-'};
for q = 1:4
    subplot(2,2,pos(q))
    hold on
    if q == 1
        M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
    elseif q == 2
        M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
    elseif q == 3
        M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
    elseif q == 4
        M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
    end
    %% arc
    R1=[[-5 5];[0,0]];
    %     R1=[[-12 12];[0,0]];
    R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',1)
    
    R1=[[0,0];[-5 5]];
    %     R1=[[0,0];[-12 12]];
    R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',1)
    
    % plot arc
    theta1 = pi/4;theta2 = 3*pi/4;the = theta1:pi/180:theta2;
    r = 10.2;x0 = 0;y0 = -10.2;x = x0 + r*cos(the);
    y = y0 + r*sin(the);
    R1=[x;y];
    R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',1)
    for c = 1:4 
        X = J_tran_group{q,c}(:,24);
        Y = J_tran_group{q,c}(:,25);
        
        %         R1=[X,Y];
        %         R2=M*R1';%旋转后坐标
        %         scatter(R2(1,:),R2(2,:),10,[0.06275	0.30588	0.5451], 'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.5)
        
        [~,~,~,X0,Y0,r_ellipse] = ErrorEllipse([X Y],0.95,colo,'-','k',1,0);
        
        R1=[r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0];
        R2=M*R1';%旋转后坐标
        plot(R2(1,:),R2(2,:),'Color',colo(c,:),'LineStyle',linestyle{c},'LineWidth',1.5)
        
        R1=[X0;Y0];
        R2=M*R1;%旋转后坐标
        scatter(R2(1),R2(2),30,'b','d','MarkerEdgeColor',colo(c,:),'LineWidth',1.5)
        
        %%
%             axis([-6.5 6.5 -6.5 6.5])
        axis([-7 7 -7 7])
%         set(gca, 'xtick',[-6 -3 0 3 6],'ytick',[-6 -3 0 3 6],'FontSize',15)
        axis square
        box off
        axis off;
        set(gca,'LineWidth',2)
%         xlabel('Tangential error (visual deg)','FontSize',15);
%         ylabel('Radial error (visual deg)','FontSize',15);
        
%         if q == 3
%             plot([-6 -1],[-6 -6],'k','LineWidth',1)
%             text(-6.5,-8,'5 cm','FontSize',10)
%         end
    end
end

print(1, '-dpng', '-r600', [datestr(now,29) '_' filename '_touchpoint_rotate'])
close all
