%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
clear all;clc;close all;
%% dMRI parameter settings
save_path='/home/user/results';
workspace_file = fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);
demographic_file='/home/user/subject_table.xlsx';
[num, txt, raw] = xlsread(demographic_file);
age=[raw{2:end,4}]';
addage=1; % 1 - Use age as a potential variable in step-wise linear regression; 0 - Not use age 
useyoung=0; % 1 - use young reproducibility group in the analysis; 0 - Not use them
%% Indexes and lsit of variables significant between-group differences for 60 dMRI acquisitions used in Labounek et al. (2020) Scientific Reports
var_indxs = [3 49 23 43 51 12 11 13 33 44 16 15 46 36 35 20 19 40 39 70 71 72 73 74 75];
table_names = {'FAwS' 'FAwK' 'f1wS' 'f1wSK' 'f1wK' 'MDwM' 'MDwm' 'MDwS' 'dwS' 'dwSK' ...
        'MDgM' 'MDgm' 'MDgSK' 'dgM' 'dgm' 'MDwgM' 'MDwgm' 'dwgM' 'dwgm' ...
        'FAwH' 'f1wH' 'MDwH' 'dwH' 'MDgH' 'dgH'};
%% HARDI-ZOOMit Interp stepwise linear regression
sbj_orig = sbj;
Y = zeros(size(table_ZOOMit_Int,1),1);
Y( [subject{:,2}]==2 | [subject{:,2}]==4, 1) = 0.5;
Y( [subject{:,2}]==1 | [subject{:,2}]==3, 1) = -0.5;

table_test = table_ZOOMit_Int(:,var_indxs);
if addage==1
        table_test(:,end+1) = age;
        table_names{1,end+1}='AGE';
end
if useyoung == 0
        table_test(sbj==3,:)=[];
        Y([subject{:,2}]==3) = [];
        sbj(sbj==3)=[];
end
table_test_norm = ( table_test - repmat(mean(table_test),size(table_test,1),1) ) ./  repmat(std(table_test),size(table_test,1),1);

X = table_test_norm;
% X(:,19) = randn(size(X,1),1);
XX = X;
Xout = [];
crl = 0;
Fmax = 0;
step = 1;
while step<=1
        [temp.b,temp.se,temp.pval,temp.inmodel,temp.stats,temp.nextstep,temp.history] =stepwisefit(XX,Y,'penter',0.059);
        temp.Y_predict = sum( repmat(temp.b(temp.inmodel==1)',size(X,1),1).*X(:,temp.inmodel==1) ,2);
        temp.crl = corrcoef([Y temp.Y_predict]);temp.crl = temp.crl(1,2);
        if abs(temp.crl) > abs(crl) && temp.stats.fstat > Fmax
                crl = temp.crl;
                Fmax =  temp.stats.fstat;
                Y_predict=temp.Y_predict;
                b = temp.b;
                se = temp.se;
                pval = temp.pval;
                inmodel = temp.inmodel;
                stats = temp.stats;
                nextstep = temp.nextstep;
                history = temp.history;
        end
        temp.pval(temp.inmodel==0) = -5;
        X1stout=find(temp.pval==max(temp.pval));
        Xout = [Xout; X1stout];
        XX(:,X1stout) = randn(size(XX,1),1)*randn(1,1);
        step = step + 1;
end
% [b,se,pval,inmodel,stats,nextstep,history] =stepwisefit(X,Y);
% Y_predict = sum( repmat(b(inmodel==1)',size(X,1),1).*X(:,inmodel==1) ,2);
% crl = corrcoef([Y Y_predict]);crl = crl(1,2);
pvl = pval(inmodel==1);beta = b(inmodel==1);SStotal = stats.SStotal;SSresid = stats.SSresid;Fvl = stats.fstat;Pvl = stats.pval;rmse=stats.rmse;
inxps = find(inmodel==1);
R2 = 1 - SSresid / SStotal;

%% HARDI-ZOOMit Interp k-means clustering
cls_results = kmeans(table_test_norm(:,inxps),2);
% cls_results = kmeans(table_test_norm(:,:),2);
cls_count = 1;
cls_num = 1;
for inter = 2:5000
        tmpidx =  kmeans(table_test_norm(:,inxps),2);
%         tmpidx =  kmeans(table_test_norm(:,:),2);
        tmpidx2 = zeros(size(tmpidx));
        tmpidx2(tmpidx==1)=2;
        tmpidx2(tmpidx==2)=1;
        cls_same = find(sum(abs(cls_results-tmpidx))==0);
        cls_same2 = find(sum(abs(cls_results-tmpidx2))==0);
        if isempty(cls_same) && ~isempty(cls_same2)
                cls_same = cls_same2;
                tmpidx = tmpidx2;
        end
        if isempty(cls_same) && isempty(cls_same2)
                cls_num = cls_num +1;
                cls_results(:,cls_num) = tmpidx;
                cls_count(1,cls_num) = 1;
        else
                cls_count(1,cls_same) = cls_count(1,cls_same) + 1;
        end
end
cls_count(cls_count==max(cls_count)) = 1;
cls_best = find(cls_count==max(cls_count));
clustered_indexes = cls_results(:,cls_best);
if clustered_indexes(1,1) == 1
        clustered_indexes(clustered_indexes==1) = 3;
        clustered_indexes(clustered_indexes==2) =1;
        clustered_indexes(clustered_indexes==3) = 2;
end

%% HARDI-ZOOMit Interp sensitivity and specificity evaluation
sbj_cls = zeros(size(Y,1),1);
for ind = 1:size(Y,1)
        if ( sbj(ind,1) == 1  || sbj(ind,1) == 3 ) && clustered_indexes(ind,1) == 1
                sbj_cls(ind,1) = sbj(ind,1);
        elseif sbj(ind,1) == 1 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 5;
        elseif sbj(ind,1) == 3 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 7;
        elseif ( sbj(ind,1) == 2  || sbj(ind,1) == 4 ) && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = sbj(ind,1);
        elseif sbj(ind,1) == 2 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 6;
        elseif sbj(ind,1) == 4 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 8;
        end
end

TP = sum(sbj_cls==2 | sbj_cls==4);
FN = sum(sbj_cls==6 | sbj_cls==8);
TN = sum(sbj_cls==1 | sbj_cls==3);
FP = sum(sbj_cls==5 | sbj_cls==7);
sensitivity = TP / (TP+FN);
specificity = TN/(TN+FP);
%% HARDI-ZOOMit Interp result visualization
% figure;plot(table_test_norm);hold on;plot(mean(table_test_norm'),'k','LineWidth',2);plot(Y,'r','LineWidth',2);plot(Y_predict,'b','LineWidth',2);hold off
HHH(1).fig = figure(2);set(HHH(1).fig,'Position',[50 50 1650 1200])
subplot(3,3,[1 2]);plot(table_test_norm,'Color',[0 1 1]);hold on;plot(Y,'r','LineWidth',4);plot(Y_predict,'b','LineWidth',3);hold off;axis([1 size(X,1) -5.3 5.3])
text(1.5,4.5,'Best model fit: ','FontWeight','bold','FontSize',12)
text(10,4.3,['Y = \beta_0 ' num2str(beta(1,1),'%+5.3f') '*' table_names{1,inxps(1,1)} ' ' ...
        num2str(beta(2,1),'%+5.3f') '*' table_names{1,inxps(1,2)} ' ' ...
        '+ \epsilon = \beta_0 + Y_p + \epsilon' ],'FontSize',12)
text(2,3.35,['Y-Y_p Pearson correlation: ' num2str(crl,'%5.3f') ],'FontSize',12)
text(11,-2.25,'\beta p-values: ','FontWeight','bold','FontSize',12)
text(17,-2.4,['\beta_1' num2str(pvl(1,1),'%.2e') '; \beta_2 ' num2str(pvl(2,1),'%.2e')],'FontSize',12)
text(11,-3.25,'Model stats: ','FontWeight','bold','FontSize',12)
text(17.2,-3.25,['F-val ' num2str(Fvl,'%5.2f') '; p-val ' num2str(Pvl,'%.2e') '; RMSE ' num2str(rmse,'%5.3f') ],'FontSize',12)
text(17.5,-4.40,[' R^2 = 1 - (SS_{\epsilon}/SS_{tot}) = ' num2str(R2*100,'%5.2f') '%' ],'FontSize',12)
title('HARDI-ZOOMit Interp','FontSize',14)
xlabel('Subject number')
ylabel('Normalized parameter')
set(gca,'FontSize',14,'LineWidth',2)

scatter_data = table_test;
MarSize=70;
subplot(3,3,3)
if sum(inmodel)==2
        if useyoung==1
                scatter( scatter_data(sbj==3 & sbj_cls==3,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,2)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', 'o')
                hold on
        end
        if sum(sbj==3 & sbj_cls==7) ~= 0
                scatter( scatter_data(sbj==3 & sbj_cls==7,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,2)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', '^')
        end
        scatter( scatter_data(sbj==1 & sbj_cls==1,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,2)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', 'o')
        if useyoung==0
                hold on
        end
        if sum(sbj==1 & sbj_cls==5) ~= 0
                scatter( scatter_data(sbj==1 & sbj_cls==5,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,2)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', '^')
        end
        scatter( scatter_data(sbj==2 & sbj_cls==2,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,2)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', '^')
        if sum(sbj==2 & sbj_cls==6) ~= 0
                scatter( scatter_data(sbj==2 & sbj_cls==6,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,2)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', 'o')
        end
        scatter( scatter_data(sbj==4 & sbj_cls==4,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,2)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', '^')
        if sum(sbj==4 & sbj_cls==8) ~= 0
                scatter( scatter_data(sbj==4 & sbj_cls==8,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,2)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', 'o')
        end
else
        if useyoung==1
                scatter3( scatter_data(sbj==3 & sbj_cls==3,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,2)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', 'o')
                hold on
        end
        if sum(sbj==3 & sbj_cls==7) ~= 0
                scatter3( scatter_data(sbj==3 & sbj_cls==7,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,2)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', '^')
        end
        scatter3( scatter_data(sbj==1 & sbj_cls==1,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,2)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', 'o')
        if useyoung==0
                hold on
        end
        if sum(sbj==1 & sbj_cls==5) ~= 0
                scatter3( scatter_data(sbj==1 & sbj_cls==5,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,2)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', '^')
        end
        scatter3( scatter_data(sbj==2 & sbj_cls==2,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,2)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', '^')
        if sum(sbj==2 & sbj_cls==6) ~= 0
                scatter3( scatter_data(sbj==2 & sbj_cls==6,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,2)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', 'o')
        end
        scatter3( scatter_data(sbj==4 & sbj_cls==4,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,2)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', '^')
        if sum(sbj==4 & sbj_cls==8) ~= 0
                scatter3( scatter_data(sbj==4 & sbj_cls==8,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,2)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', 'o')
        end
end
hold off
xlabel(table_names{1,inxps(1,1)})
ylabel(table_names{1,inxps(1,2)})
if sum(inmodel)>2
        zlabel(table_names{1,inxps(1,3)})
        az = 57.76;
        el = 28.56;
        view(az, el);
else
        grid on
end
% title({'Data projections in spaces';'of significant variables'},'FontSize',14)
title({'K-means';['sensitivity = ' num2str(sensitivity*100,'%5.2f') '%; specificity = ' num2str(specificity*100,'%5.2f') '%']},'FontSize',12)
set(gca,'FontSize',12,'LineWidth',2)

%% DTI-RESOLVE Non-Interp stepwise linear regression
table_test = table_RESOLVE(:,var_indxs);
if addage==1
        table_test(:,end+1) = age;
end
sbj=sbj_orig;
if useyoung == 0
        table_test(sbj==3,:)=[];
        sbj(sbj==3)=[];
end
table_test_norm = ( table_test - repmat(mean(table_test),size(table_test,1),1) ) ./  repmat(std(table_test),size(table_test,1),1);

X = table_test_norm;
XX = X;
Xout = [];
crl = 0;
Fmax = 0;
step = 1;
while step<=1
        [temp.b,temp.se,temp.pval,temp.inmodel,temp.stats,temp.nextstep,temp.history] =stepwisefit(XX,Y,'penter',0.052);
        temp.Y_predict = sum( repmat(temp.b(temp.inmodel==1)',size(X,1),1).*X(:,temp.inmodel==1) ,2);
        temp.crl = corrcoef([Y temp.Y_predict]);temp.crl = temp.crl(1,2);
        if abs(temp.crl) > abs(crl) && temp.stats.fstat > Fmax
                crl = temp.crl;
                Fmax =  temp.stats.fstat;
                Y_predict=temp.Y_predict;
                b = temp.b;
                se = temp.se;
                pval = temp.pval;
                inmodel = temp.inmodel;
                stats = temp.stats;
                nextstep = temp.nextstep;
                history = temp.history;
        end
        temp.pval(temp.inmodel==0) = -5;
        X1stout=find(temp.pval==max(temp.pval));
        Xout = [Xout; X1stout];
        XX(:,X1stout) = randn(size(XX,1),1)*randn(1,1);
        step = step + 1;
end
% [b,se,pval,inmodel,stats,nextstep,history] =stepwisefit(X,Y);
% Y_predict = sum( repmat(b(inmodel==1)',size(X,1),1).*X(:,inmodel==1) ,2);
% crl = corrcoef([Y Y_predict]);crl = crl(1,2);
pvl = pval(inmodel==1);beta = b(inmodel==1);SStotal = stats.SStotal;SSresid = stats.SSresid;Fvl = stats.fstat;Pvl = stats.pval;rmse=stats.rmse;
inxps = find(inmodel==1);
R2 = 1 - SSresid / SStotal;

%% DTI-RESOLVE Non-Interp k-means clustering
cls_results = kmeans(table_test_norm(:,inxps),2);
cls_count = 1;
cls_num = 1;
for inter = 2:5000
        tmpidx =  kmeans(table_test_norm(:,inxps),2);
        tmpidx2 = zeros(size(tmpidx));
        tmpidx2(tmpidx==1)=2;
        tmpidx2(tmpidx==2)=1;
        cls_same = find(sum(abs(cls_results-tmpidx))==0);
        cls_same2 = find(sum(abs(cls_results-tmpidx2))==0);
        if isempty(cls_same) && ~isempty(cls_same2)
                cls_same = cls_same2;
                tmpidx = tmpidx2;
        end
        if isempty(cls_same) && isempty(cls_same2)
                cls_num = cls_num +1;
                cls_results(:,cls_num) = tmpidx;
                cls_count(1,cls_num) = 1;
        else
                cls_count(1,cls_same) = cls_count(1,cls_same) + 1;
        end
end
% cls_count(cls_count==max(cls_count)) = 1;
cls_best = find(cls_count==max(cls_count));
clustered_indexes = cls_results(:,cls_best);
if clustered_indexes(1,1) == 1
        clustered_indexes(clustered_indexes==1) = 3;
        clustered_indexes(clustered_indexes==2) =1;
        clustered_indexes(clustered_indexes==3) = 2;
end

%% DTI-RESOLVE Non-Interp sensitivity and specificity visualization
sbj_cls = zeros(size(sbj,1),1);
for ind = 1:size(sbj,1)
        if ( sbj(ind,1) == 1  || sbj(ind,1) == 3 ) && clustered_indexes(ind,1) == 1
                sbj_cls(ind,1) = sbj(ind,1);
        elseif sbj(ind,1) == 1 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 5;
        elseif sbj(ind,1) == 3 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 7;
        elseif ( sbj(ind,1) == 2  || sbj(ind,1) == 4 ) && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = sbj(ind,1);
        elseif sbj(ind,1) == 2 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 6;
        elseif sbj(ind,1) == 4 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 8;
        end
end

TP = sum(sbj_cls==2 | sbj_cls==4);
FN = sum(sbj_cls==6 | sbj_cls==8);
TN = sum(sbj_cls==1 | sbj_cls==3);
FP = sum(sbj_cls==5 | sbj_cls==7);
sensitivity = TP / (TP+FN);
specificity = TN/(TN+FP);

%% DTI-RESOLVE Non-Interp result visualization
% figure;plot(table_test_norm);hold on;plot(mean(table_test_norm'),'k','LineWidth',2);plot(Y,'r','LineWidth',2);plot(Y_predict,'b','LineWidth',2);hold off
subplot(3,3,[4 5]);plot(table_test_norm,'Color',[0 1 1]);hold on;plot(Y,'r','LineWidth',4);plot(Y_predict,'b','LineWidth',3);hold off;axis([1 size(X,1) -5.3 5.3])
text(1.5,4.5,'Best model fit: ','FontWeight','bold','FontSize',12)
if sum(inmodel)==2
        text(10,4.3,['Y = \beta_0 ' num2str(beta(1,1),'%+5.3f') '*' table_names{1,inxps(1,1)} ' ' ...
                num2str(beta(2,1),'%+5.3f') '*' table_names{1,inxps(1,2)} ' ' ...
                '+ \epsilon = \beta_0 + Y_p + \epsilon' ],'FontSize',12)
else
        text(10,4.3,['Y = \beta_0 ' num2str(beta(1,1),'%+5.3f') '*' table_names{1,inxps(1,1)} ' ' ...
                num2str(beta(2,1),'%+5.3f') '*' table_names{1,inxps(1,2)} ' ' ...
                num2str(beta(3,1),'%+5.3f') '*' table_names{1,inxps(1,3)} ' ' ...
                '+ \epsilon = \beta_0 + Y_p + \epsilon' ],'FontSize',12)
end
text(2,3.35,['Y-Y_p Pearson correlation: ' num2str(crl,'%5.3f') ],'FontSize',12)
text(11,-2.25,'\beta p-values: ','FontWeight','bold','FontSize',12)
text(17,-2.4,['\beta_1' num2str(pvl(1,1),'%.2e') '; \beta_2 ' num2str(pvl(2,1),'%.2e')  '; \beta_3 ' num2str(pvl(3,1),'%.2e') ],'FontSize',12)
text(11,-3.25,'Model stats: ','FontWeight','bold','FontSize',12)
text(17.2,-3.25,['F-val ' num2str(Fvl,'%5.2f') '; p-val ' num2str(Pvl,'%.2e') '; RMSE ' num2str(rmse,'%5.3f') ],'FontSize',12)
text(17.5,-4.40,[' R^2 = 1 - (SS_{\epsilon}/SS_{tot}) = ' num2str(R2*100,'%5.2f') '%' ],'FontSize',12)
title('DTI-RESOLVE','FontSize',14)
xlabel('Subject number')
ylabel('Normalized parameter')
set(gca,'FontSize',14,'LineWidth',2)


scatter_data = table_test;
subplot(3,3,6)
if sum(inmodel)==2
        if useyoung==1
                scatter( scatter_data(sbj==3 & sbj_cls==3,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,2)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', 'o')
                hold on
        end
        if sum(sbj==3 & sbj_cls==7) ~= 0
                scatter( scatter_data(sbj==3 & sbj_cls==7,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,2)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', '^')
        end
        scatter( scatter_data(sbj==1 & sbj_cls==1,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,2)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', 'o')
        if useyoung==0
                hold on
        end
        if sum(sbj==1 & sbj_cls==5) ~= 0
                scatter( scatter_data(sbj==1 & sbj_cls==5,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,2)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', '^')
        end
        scatter( scatter_data(sbj==2 & sbj_cls==2,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,2)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', '^')
        if sum(sbj==2 & sbj_cls==6) ~= 0
                scatter( scatter_data(sbj==2 & sbj_cls==6,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,2)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', 'o')
        end
        scatter( scatter_data(sbj==4 & sbj_cls==4,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,2)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', '^')
        if sum(sbj==4 & sbj_cls==8) ~= 0
                scatter( scatter_data(sbj==4 & sbj_cls==8,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,2)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', 'o')
        end
else
        if useyoung==1
                scatter3( scatter_data(sbj==3 & sbj_cls==3,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,2)), scatter_data(sbj==3 & sbj_cls==3,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', 'o')
                hold on
        end
        if sum(sbj==3 & sbj_cls==7) ~= 0
                scatter3( scatter_data(sbj==3 & sbj_cls==7,inxps(1,1)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,2)), scatter_data(sbj==3 & sbj_cls==7,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', '^')
        end
        scatter3( scatter_data(sbj==1 & sbj_cls==1,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,2)), scatter_data(sbj==1 & sbj_cls==1,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', 'o')
        if useyoung==0
                hold on
        end
        if sum(sbj==1 & sbj_cls==5) ~= 0
                scatter3( scatter_data(sbj==1 & sbj_cls==5,inxps(1,1)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,2)), scatter_data(sbj==1 & sbj_cls==5,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', '^')
        end
        scatter3( scatter_data(sbj==2 & sbj_cls==2,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,2)), scatter_data(sbj==2 & sbj_cls==2,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', '^')
        if sum(sbj==2 & sbj_cls==6) ~= 0
                scatter3( scatter_data(sbj==2 & sbj_cls==6,inxps(1,1)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,2)), scatter_data(sbj==2 & sbj_cls==6,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', 'o')
        end
        scatter3( scatter_data(sbj==4 & sbj_cls==4,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,2)), scatter_data(sbj==4 & sbj_cls==4,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', '^')
        if sum(sbj==4 & sbj_cls==8) ~= 0
                scatter3( scatter_data(sbj==4 & sbj_cls==8,inxps(1,1)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,2)), scatter_data(sbj==4 & sbj_cls==8,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', 'o')
        end
end
hold off
xlabel(table_names{1,inxps(1,1)})
ylabel(table_names{1,inxps(1,2)})
if sum(inmodel)>2
        zlabel(table_names{1,inxps(1,3)})
        az = 57.76;
        el = 28.56;
        view(az, el);
else
        grid on
end
title({'K-means';['sensitivity = ' num2str(sensitivity*100,'%5.2f') '%; specificity = ' num2str(specificity*100,'%5.2f') '%'];' '},'FontSize',12)
set(gca,'FontSize',12,'LineWidth',2)

%% HARDI-ZOOMit Non-Interp stepwise linear regression
subject_zommnonint = [subject{:,2}];
subject_zommnonint(MD_gauss_wm_mean(:,2) == 0) = [];
sbj_NtInt = subject_zommnonint';
Y = zeros(size(subject_zommnonint,2),1);
Y( subject_zommnonint ==2 | subject_zommnonint==4, 1) = 0.5;
Y( subject_zommnonint==1 | subject_zommnonint==3, 1) = -0.5;

table_test = table_ZOOMit_NotInt(:,var_indxs);
table_test(sum(table_test,2)==0,:) = [];
if addage==1
        ag = age;
        ag(MD_gauss_wm_mean(:,2) == 0) = [];
        table_test(:,end+1) = ag;
end
if useyoung == 0
        table_test(sbj_NtInt==3,:)=[];
        Y(sbj_NtInt==3)=[];
        sbj_NtInt(sbj_NtInt==3)=[];
end
table_test_norm = ( table_test - repmat(mean(table_test),size(table_test,1),1) ) ./  repmat(std(table_test),size(table_test,1),1);

X = table_test_norm;
XX = X;
Xout = [];
crl = 0;
Fmax = 0;
step = 1;
while step<=1
        [temp.b,temp.se,temp.pval,temp.inmodel,temp.stats,temp.nextstep,temp.history] =stepwisefit(XX,Y,'penter',0.12);
        temp.Y_predict = sum( repmat(temp.b(temp.inmodel==1)',size(X,1),1).*X(:,temp.inmodel==1) ,2);
        temp.crl = corrcoef([Y temp.Y_predict]);temp.crl = temp.crl(1,2);
        if abs(temp.crl) > abs(crl) && temp.stats.fstat > Fmax
                crl = temp.crl;
                Fmax =  temp.stats.fstat;
                Y_predict=temp.Y_predict;
                b = temp.b;
                se = temp.se;
                pval = temp.pval;
                inmodel = temp.inmodel;
                stats = temp.stats;
                nextstep = temp.nextstep;
                history = temp.history;
        end
        temp.pval(temp.inmodel==0) = -5;
        X1stout=find(temp.pval==max(temp.pval));
        Xout = [Xout; X1stout];
        XX(:,X1stout) = randn(size(XX,1),1)*randn(1,1);
        step = step + 1;
end
% [b,se,pval,inmodel,stats,nextstep,history] =stepwisefit(X,Y);
% Y_predict = sum( repmat(b(inmodel==1)',size(X,1),1).*X(:,inmodel==1) ,2);
% crl = corrcoef([Y Y_predict]);crl = crl(1,2);
pvl = pval(inmodel==1);beta = b(inmodel==1);SStotal = stats.SStotal;SSresid = stats.SSresid;Fvl = stats.fstat;Pvl = stats.pval;rmse=stats.rmse;
inxps = find(inmodel==1);
R2 = 1 - SSresid / SStotal;

%% HARDI-ZOOMit Non-Interp k-means clustering
cls_results = kmeans(table_test_norm(:,inxps),2);
cls_count = 1;
cls_num = 1;
for inter = 2:5000
        tmpidx =  kmeans(table_test_norm(:,inxps),2);
        tmpidx2 = zeros(size(tmpidx));
        tmpidx2(tmpidx==1)=2;
        tmpidx2(tmpidx==2)=1;
        cls_same = find(sum(abs(cls_results-tmpidx))==0);
        cls_same2 = find(sum(abs(cls_results-tmpidx2))==0);
        if isempty(cls_same) && ~isempty(cls_same2)
                cls_same = cls_same2;
                tmpidx = tmpidx2;
        end
        if isempty(cls_same) && isempty(cls_same2)
                cls_num = cls_num +1;
                cls_results(:,cls_num) = tmpidx;
                cls_count(1,cls_num) = 1;
        else
                cls_count(1,cls_same) = cls_count(1,cls_same) + 1;
        end
end
% cls_count(cls_count==max(cls_count)) = 1;
% cls_count(cls_count==max(cls_count))=0; %!!!!!!! Taking 2nd best kmenas result
cls_best = find(cls_count==max(cls_count));
clustered_indexes = cls_results(:,cls_best);
if clustered_indexes(1,1) == 1
        clustered_indexes(clustered_indexes==1) = 3;
        clustered_indexes(clustered_indexes==2) =1;
        clustered_indexes(clustered_indexes==3) = 2;
end

%% HARDI-ZOOMit Non-Interp sensitivity and specificity evaluation
sbj_cls = zeros(size(sbj_NtInt,1),1);
for ind = 1:size(sbj_NtInt,1)
        if ( sbj_NtInt(ind,1) == 1  || sbj_NtInt(ind,1) == 3 ) && clustered_indexes(ind,1) == 1
                sbj_cls(ind,1) = sbj_NtInt(ind,1);
        elseif sbj_NtInt(ind,1) == 1 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 5;
        elseif sbj_NtInt(ind,1) == 3 && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = 7;
        elseif ( sbj_NtInt(ind,1) == 2  || sbj_NtInt(ind,1) == 4 ) && clustered_indexes(ind,1) == 2
                sbj_cls(ind,1) = sbj_NtInt(ind,1);
        elseif sbj_NtInt(ind,1) == 2 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 6;
        elseif sbj_NtInt(ind,1) == 4 && clustered_indexes(ind,1) ==1
                sbj_cls(ind,1) = 8;
        end
end

TP = sum(sbj_cls==2 | sbj_cls==4);
FN = sum(sbj_cls==6 | sbj_cls==8);
TN = sum(sbj_cls==1 | sbj_cls==3);
FP = sum(sbj_cls==5 | sbj_cls==7);
sensitivity = TP / (TP+FN);
specificity = TN/(TN+FP);

%% HARDI-ZOOMit Non-Interp result visualization
subplot(3,3,[7 8]);LL1=plot(table_test_norm(:,1),'Color',[0 1 1]);hold on;plot(table_test_norm(:,2:end),'Color',[0 1 1]);LL2=plot(Y,'r','LineWidth',4);LL3=plot(Y_predict,'b','LineWidth',3);hold off;axis([1 46 -5.3 5.3])
text(2,4.5,'Best model fit: ','FontWeight','bold','FontSize',12)
text(10,4.3,['Y = \beta_0 ' num2str(beta(1,1),'%+5.3f') '*' table_names{1,inxps(1,1)} ' ' ...
        num2str(beta(2,1),'%+5.3f') '*' table_names{1,inxps(1,2)} ' ' ...
        num2str(beta(3,1),'%+5.3f') '*' table_names{1,inxps(1,3)} ' ' ...
        num2str(beta(4,1),'%+5.3f') '*' table_names{1,inxps(1,4)} '+ \epsilon = \beta_0 + Y_p + \epsilon' ],'FontSize',12)
text(2,3.35,['Y-Y_p Pearson correlation: ' num2str(crl,'%5.3f') ],'FontSize',12)
text(10,-2.25,'\beta p-values: ','FontWeight','bold','FontSize',12)
text(17,-2.4,['\beta_1' num2str(pvl(1,1),'%.2e') '; \beta_2 ' num2str(pvl(2,1),'%.2e') '; \beta_3 ' num2str(pvl(3,1),'%.2e')   '; \beta_4 ' num2str(pvl(4,1),'%.2e')],'FontSize',12)
text(10,-3.25,'Model stats: ','FontWeight','bold','FontSize',12)
text(17.2,-3.25,['F-val ' num2str(Fvl,'%5.2f') '; p-val ' num2str(Pvl,'%.2e') '; RMSE ' num2str(rmse,'%5.3f') ],'FontSize',12)
text(17.5,-4.40,[' R^2 = 1 - (SS_{\epsilon}/SS_{tot}) = ' num2str(R2*100,'%5.2f') '%' ],'FontSize',12)
title('HARDI-ZOOMit Non-Interp','FontSize',14)
xlabel('Subject number')
ylabel('Normalized parameter')
legend([LL1 LL2 LL3],{'Normalized dMRI-derived parameters N(0,1)','Y','Y_P'})
set(gca,'FontSize',14,'LineWidth',2)

% inxps2=inxps;
% inxps=inxps(2:end);
% % inxps(pvl==max(pvl))=[];
scatter_data = table_test;
subplot(3,3,9)
if useyoung==1
        scatter3( scatter_data(sbj_NtInt==3 & sbj_cls==3,inxps(1,1)), scatter_data(sbj_NtInt==3 & sbj_cls==3,inxps(1,4)), scatter_data(sbj_NtInt==3 & sbj_cls==3,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', 'o')
        hold on
end
if sum(sbj_NtInt==3 & sbj_cls==7) ~= 0
        scatter3( scatter_data(sbj_NtInt==3 & sbj_cls==7,inxps(1,1)), scatter_data(sbj_NtInt==3 & sbj_cls==7,inxps(1,4)), scatter_data(sbj_NtInt==3 & sbj_cls==7,inxps(1,3)),MarSize,'MarkerEdgeColor','k','LineWidth',2,'Marker', '^')
end
scatter3( scatter_data(sbj_NtInt==1 & sbj_cls==1,inxps(1,1)), scatter_data(sbj_NtInt==1 & sbj_cls==1,inxps(1,4)), scatter_data(sbj_NtInt==1 & sbj_cls==1,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', 'o')
if useyoung==0
        hold on
end
if sum(sbj_NtInt==1 & sbj_cls==5) ~= 0
        scatter3( scatter_data(sbj_NtInt==1 & sbj_cls==5,inxps(1,1)), scatter_data(sbj_NtInt==1 & sbj_cls==5,inxps(1,4)), scatter_data(sbj_NtInt==1 & sbj_cls==5,inxps(1,3)),MarSize,'MarkerEdgeColor','r','LineWidth',2,'Marker', '^')
end
scatter3( scatter_data(sbj_NtInt==2 & sbj_cls==2,inxps(1,1)), scatter_data(sbj_NtInt==2 & sbj_cls==2,inxps(1,4)), scatter_data(sbj_NtInt==2 & sbj_cls==2,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', '^')
if sum(sbj_NtInt==2 & sbj_cls==6) ~= 0
        scatter3( scatter_data(sbj_NtInt==2 & sbj_cls==6,inxps(1,1)), scatter_data(sbj_NtInt==2 & sbj_cls==6,inxps(1,4)), scatter_data(sbj_NtInt==2 & sbj_cls==6,inxps(1,3)),MarSize,'MarkerEdgeColor',[0 0.749 1],'LineWidth',2,'Marker', 'o')
end
scatter3( scatter_data(sbj_NtInt==4 & sbj_cls==4,inxps(1,1)), scatter_data(sbj_NtInt==4 & sbj_cls==4,inxps(1,4)), scatter_data(sbj_NtInt==4 & sbj_cls==4,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', '^')
if sum(sbj_NtInt==4 & sbj_cls==8) ~= 0
        scatter3( scatter_data(sbj_NtInt==4 & sbj_cls==8,inxps(1,1)), scatter_data(sbj_NtInt==4 & sbj_cls==8,inxps(1,4)), scatter_data(sbj_NtInt==4 & sbj_cls==8,inxps(1,3)),MarSize,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',2,'Marker', 'o')
end
hold off
xlabel(table_names{1,inxps(1,1)})
ylabel(table_names{1,inxps(1,4)})
zlabel(table_names{1,inxps(1,3)})
az = 57.76;
el = 28.56;
view(az, el);
% title({'Data projections in spaces';'of significant variables'},'FontSize',14)
title({'K-means';['sensitivity = ' num2str(sensitivity*100,'%5.2f') '%; specificity = ' num2str(specificity*100,'%5.2f') '%']},'FontSize',12)
set(gca,'FontSize',12,'LineWidth',2)