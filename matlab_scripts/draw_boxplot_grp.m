function draw_boxplot_grp(x,group,ymin,ymax,ind,will_pvals,WillThr,triangle_pos,dot_pos)
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com

positions = [1 1.05 1.1 1.15];
positions2 = 0.15+0.15+positions;
positions3 = 0.25+0.35+positions;
positions = [positions positions2 positions3];
color = [0.65 0.65 0.65;
        0 0.749 1;
        1 0 0;
        0 0 0];
color = repmat(color,3,1);
color_reverse = [0 0 0;
        1 0 0
        0 0.749 1;
        0.65 0.65 0.65];
color_reverse = repmat(color_reverse,3,1);

bh = boxplot(x,group, 'positions', positions);

pos1=mean(positions(1:4));
pos2=mean(positions(5:8));
pos3=mean(positions(9:12));
set(gca,'xtick',[pos1 pos2 pos3])
set(gca,'xticklabel',{'HZi','HZni','DR'})

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.65);
end
for i=1:size(bh,1) % <- # graphics handles/x
        for j = 1:size(bh,2)
               if i == 1 || i == 2
                        set(bh(i,j),'linewidth',2,'Color',color_reverse(j,:),'LineStyle','-');
%                         set(bh(i,j),'Color',color_reverse(j,:));
               elseif i == 3 || i == 4
                       set(bh(i,j),'linewidth',3,'Color','k','LineStyle','-');
               elseif i == 7
                       set(bh(i,j),'Marker','.','MarkerSize',9);
                elseif i ~= 6
                        set(bh(i,j),'linewidth',2,'Color',color_reverse(j,:));
                else
                        set(bh(i,j),'linewidth',4,'Color',[0 0 0]);
                end
%                 disp([num2str(i) '    ' num2str(j)])
                %                 disp(sprintf('working on component: %3d = %s',i,get(bh(i,1),'tag')));
                pause(.05);
        end
end
hold on
if will_pvals(1,8) < WillThr
        plot(positions(4),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(1,:))
end
if will_pvals(1,7) < WillThr
        plot(positions(3),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(2,:))
end
if will_pvals(1,10) < WillThr
        plot(positions(8),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(1,:))
end
if will_pvals(1,9) < WillThr
        plot(positions(7),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(2,:))
end
if will_pvals(1,12) < WillThr
        plot(positions(12),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(1,:))
end
if will_pvals(1,11) < WillThr
        plot(positions(11),triangle_pos,'^','LineStyle','none','LineWidth',2,'MarkerSize',7,'Color',color(2,:))
end
if will_pvals(1,2) < WillThr
        plot(positions(4),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(1,:))
end
if will_pvals(1,1) < WillThr
        plot(positions(3),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(2,:))
end
if will_pvals(1,4) < WillThr
        plot(positions(8),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(1,:))
end
if will_pvals(1,3) < WillThr
        plot(positions(7),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(2,:))
end
if will_pvals(1,6) < WillThr
        plot(positions(12),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(1,:))
end
if will_pvals(1,5) < WillThr
        plot(positions(11),dot_pos,'*','LineStyle','none','LineWidth',3,'MarkerSize',5,'Color',color(2,:))
end
if ind == 21
        plot([0 0.2],[-5 -5],'Marker','.','MarkerSize',9,'LineStyle','none','Color',color(3,:))
        c = get(gca, 'Children');
        hleg1 = legend(c(1:9),'Outlayers','C-M significance','C-S significance','R-M significance','R-S significance','Reproducibility (R)', 'Controls (C)', 'M compression patients', 'S compression patients' );
end
hold off
axis([0.9 1.85 ymin ymax])
grid on
set(gca,'FontSize',12,'LineWidth',2)