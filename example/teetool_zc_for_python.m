
%% select file from
currentFolder = pwd;
namelist = dir('**/*.mat');
figure
% set(gcf,'position',[30,300,500,400]);
hold on
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;...
    0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;...
    0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;...
    0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
linestyle = {'--','--','-','-','--','--','-','-','--','--','-','-','--','--','-','-'};
n = 1;s_name=[];
for s = 1:length(namelist)
    if ~isempty(strfind(namelist(s).name,'rotrj_ci')) %&& ~isempty(strfind(namelist(ses).name,'0815'))
        tic
        load(namelist(s).name,'x','y','z')
        s_name{n,1} = namelist(s).name;
        
        z = squeeze(z);
        [C,h]=contour(x,y,z);
        h.LineColor = colo(n,:);
        h.LineStyle = linestyle{n};
        h.Visible = 'off';
        npoints = C(2,1);
        fill(C(1,2:npoints+1),C(2,2:npoints+1),colo(n,:),'EdgeColor','none','FaceAlpha',0.15) % Ìî³ä

        n = n+1;
        toc
    end
end

axis off
axis square

%% plot mean

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
                R2=M*R1';%Ðý×ªºó×ø±ê
                plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
    end
end

set(gcf,'position',[30,300,400,300])
