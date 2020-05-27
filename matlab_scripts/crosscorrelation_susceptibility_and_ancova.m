%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
%
% Copyright 2016-2020 Rene Labounek (1,2,3,4), Jan Valosek (1,2) and Petr Hlustik (1,2)
%
% 1 - University Hospital Olomouc, Olomouc, Czech Republic
% 2 - Palacky University Olomouc, Olomouc, Czech Republic
% 3 - University Hospital Brno, Brno, Czech Republic 
% 4 - University of Minnesota, Minneapolis, US
%
% This file is part of sc-dmri-myelopathy available at: https://github.com/renelabounek/sc-dmri-myelopathy
%
% Please, cite sc-dmri-myelopathy as:
% Labounek R, Valosek J, Horak T, Svatkova A, Bednarik P, Vojtisek L, Horakova M, Nestrasil I,
% Lenglet C, Cohen-Adad J, Bednarik J and Hlustik P. HARDI-ZOOMit protocol improves specificity
% to microstructural changes in presymptomatic myelopathy. Scientific Reports [Revised; Under review]
%
% sc-dmri-myelopathy is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% sc-dmri-myelopathy is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with sc-dmri-myelopathy.  If not, see <https://www.gnu.org/licenses/>.
%
clc;
clear all;
close all;
%% dMRI parameter settings
save_path='/home/user/results';
demographic_file='/home/user/subject_table.xlsx';
thr_cor=0.05;
workspace_file = fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);
[num, txt, raw] = xlsread(demographic_file);
age=[raw{2:end,4}]';
cols = size(raw,2);
%% Indexes and lsit of variables significant between-group differences for 60 dMRI acquisitions used in Labounek et al. (2020) Scientific Reports
% Last 6 variables are measures of sscueptibility artifact levels from C3-C6, C3 and C5-C6 regions of interest
var_indxs = [3 49 23 43 51 12 11 13 33 44 16 15 46 36 35 20 19 40 39 70 71 72 73 74 75 60 61 62 63 64 65];
table_names = {'FAwS' 'FAwK' 'f1wS' 'f1wSK' 'f1wK' 'MDwM' 'MDwm' 'MDwS' 'dwS' 'dwSK' ...
        'MDgM' 'MDgm' 'MDgSK' 'dgM' 'dgm' 'MDwgM' 'MDwgm' 'dwgM' 'dwgm' ...
        'FAwH' 'f1wH' 'MDwH' 'dwH' 'MDgH' 'dgH' ...
        'S36M' 'S3M' 'S56M' 'S36S' 'S3S' 'S56S'};
varsignum= 25;
%% HARDI-ZOOMit Interp cross-corrrelations and ANCOVA
X = table_ZOOMit_Int(:,var_indxs);
[RHZi, p_RHZi] = corrcoef(X,'Rows','pairwise');
RHZi = triu(RHZi);
p_RHZi = triu(p_RHZi);

p_ancova_RHZI=zeros(2,varsignum);
% Confounding variable
confounder=age;
% Grouping variable
group = [subject{:,2}]';
group(group == 1 | group == 3) = 0;
group(group == 2 | group == 4) = 1;
% group vocabulary:
% 3 - reproducibility (yound healthy subjects acquired twice with inter-scan interval >1 day and <29 weeks)
% 1 - age-comparable healthy controls
% 2 - mild asymptomatic degenerative cervical cord compression patients
% 4 - severe asymptomatic degenerative cervical cord compression patients
for metid = 1:varsignum
        p_tmp = anovan(table_ZOOMit_Int(:,var_indxs(1,metid)),{group,confounder},'Continuous',2,'varnames',{'Group','confounder'},'display','off');
        p_ancova_RHZI(1,metid)= p_tmp(1);
end

group = [subject{:,2}]';
group(group == 1) = 0;
group(group == 2 | group == 4) = 1;
confounder=age(group~=3);
data = table_ZOOMit_Int(:,var_indxs(1:varsignum));
data(group==3,:) = [];
group(group==3) = [];
for metid = 1:varsignum
        p_tmp = anovan(data(:,metid),{group,confounder},'Continuous',2,'varnames',{'Group','confounder'},'display','off');
        p_ancova_RHZI(2,metid)= p_tmp(1);
end
%% DTI-RESOLVE cross-corrrelations and ANCOVA
X = table_RESOLVE(:,var_indxs);
[RDR, p_RDR]= corrcoef(X,'Rows','pairwise');
RDR = triu(RDR);
p_RDR = triu(p_RDR);

p_ancova_RDR=zeros(2,varsignum);
confounder=age;
group = [subject{:,2}]';
group(group == 1 | group == 3) = 0;
group(group == 2 | group == 4) = 1;
for metid = 1:varsignum
        p_tmp = anovan(table_RESOLVE(:,var_indxs(1,metid)),{group,confounder},'Continuous',2,'varnames',{'Group','confounder'},'display','off');
%         p_tmp = anovan(table_ZOOMit_Int(:,var_indxs(1,metid)),{group,confounder},'Continuous',2,'varnames',{'Group','confounder'});
        p_ancova_RDR(1,metid)= p_tmp(1);
end

group = [subject{:,2}]';
group(group == 1) = 0;
group(group == 2 | group == 4) = 1;
confounder=age(group~=3);
data = table_RESOLVE(:,var_indxs(1:varsignum));
data(group==3,:) = [];
group(group==3) = [];
for metid = 1:varsignum
        p_tmp = anovan(data(:,metid),{group,confounder},'Continuous',2,'varnames',{'Group','confounder'},'display','off');
        p_ancova_RDR(2,metid)= p_tmp(1);
end
%% HARDI-ZOOMit Non-Interp cross-correlations
X = table_ZOOMit_NotInt(:,var_indxs);
X(sum(X,2)==0,:) = [];
[RHZni, p_RHZni] = corrcoef(X,'Rows','pairwise');
RHZni = triu(RHZni);
p_RHZni = triu(p_RHZni);
%% mean iduced cross-correlations between dMRI and off-resonance effects
temp = abs(RHZi(1:25,26:end));
offres_crosscorr(1,1) = mean(temp(:));
offres_crosscorr(2,1) = std(temp(:));

temp = abs(RHZni(1:25,26:end));
offres_crosscorr(1,2) = mean(temp(:));
offres_crosscorr(2,2) = std(temp(:));

temp = abs(RDR(1:25,26:end));
offres_crosscorr(1,3) = mean(temp(:));
offres_crosscorr(2,3) = std(temp(:));

t_value_R(1,1) = offres_crosscorr(1,1) / sqrt( (1-offres_crosscorr(1,1)^2)/(60-2)  );
p_value_R(1,1) = 2*(1-tcdf(abs(t_value_R(1,1)),60-1));
t_value_R(1,2) = offres_crosscorr(1,2) / sqrt( (1-offres_crosscorr(1,2)^2)/(22-2)  );
p_value_R(1,2) = 2*(1-tcdf(abs(t_value_R(1,2)),22-1));
t_value_R(1,3) = offres_crosscorr(1,3) / sqrt( (1-offres_crosscorr(1,3)^2)/(60-2)  );
p_value_R(1,3) = 2*(1-tcdf(abs(t_value_R(1,3)),60-1));

offres_crosscorr = round(offres_crosscorr*100)/100;
%% Susceptibility statistics
Sarea = abs(RHZi(1:25,26:end));
mean(Sarea(:))
std(Sarea(:))
Sarea = abs(RHZni(1:25,26:end));
mean(Sarea(:))
std(Sarea(:))
Sarea = abs(RDR(1:25,26:end));
mean(Sarea(:))
std(Sarea(:))
%% Ancova result full table
% R - reproducibility
% C - age-comparable healthy controls
% M - mild asymptomatic degenerative cervical cord compression patients
% S - severe asymptomatic degenerative cervical cord compression patients
% First column are p-values of HARDI-ZOOMit Interp between R+C vs M+S
% Second column are p-values of HARDI-ZOOMit Interp between C vs M+S
% Third column are p-values of DTI-RESOLVE Non-Interp between R+C vs M+S
% Fourth column are p-values of DTI-RESOLVE Non-Interp between C vs M+S
p_ancova = [p_ancova_RHZI; p_ancova_RDR]';
%% Cross-correlation plot
[C, D] = meshgrid(1:size(RHZi,1)+1);

h(1).fig = figure(1);
set(h(1).fig, 'Position', [50, 50, 2450, 590]);

subplot(1,3,1)
temp = RHZi;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = -1;
surf(C,D,temp,'edgecolor','none')
title('HARDI-ZOOMit Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RHZi,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RHZi,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap jet
colorbar
axis([ 1 size(RHZi,1)+1 1 size(RHZi,1)+1 min(RHZi(:)) max(RHZi(:)) ])
az = 90;
el = 90;
view(az, el);

subplot(1,3,2)
temp = RDR;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = -1;
surf(C,D,temp,'edgecolor','none')
title('DTI-RESOLVE Non-Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RDR,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RDR,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap jet
colorbar
axis([ 1 size(RDR,1)+1 1 size(RDR,1)+1 min(RDR(:)) max(RDR(:)) ])
az = 90;
el = 90;
view(az, el);

subplot(1,3,3)
temp = RHZni ;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = -1;
surf(C,D,temp,'edgecolor','none')
title('HARDI-ZOOMit Non-Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RHZni ,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RHZni ,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap jet
colorbar
axis([ 1 size(RHZni ,1)+1 1 size(RHZni ,1)+1 min(RHZni (:)) max(RHZni (:)) ])
az = 90;
el = 90;
view(az, el);
%% Cross-correlation p-value plot
h(3).fig = figure(3);
% set(h(3).fig, 'Position', [50, 50, 2300, 540]);
set(h(3).fig, 'Position', [50, 50, 2450, 590]);

subplot(1,3,1)
temp = p_RHZi;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = 1;
temp(temp>0.05) = 0.05;
surf(C,D,temp,'edgecolor','none')
title('HARDI-ZOOMit Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RHZi,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RHZi,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap hot
colorbar('limits',[0 0.05])
axis([ 1 size(RHZi,1)+1 1 size(RHZi,1)+1 min(RHZi(:)) max(RHZi(:)) ])
az = 90;
el = 90;
view(az, el);

subplot(1,3,2)
temp = p_RDR;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = 1;
temp(temp>0.05) = 0.05;
surf(C,D,temp,'edgecolor','none')
title('DTI-RESOLVE Non-Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RDR,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RDR,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap hot
colorbar
axis([ 1 size(RDR,1)+1 1 size(RDR,1)+1 min(RDR(:)) max(RDR(:)) ])
az = 90;
el = 90;
view(az, el);

subplot(1,3,3)
temp = p_RHZni ;
temp(end+1,:) = 0;
temp(:,end+1) = 0;
temp(end,end) = 1;
temp(temp>0.05) = 0.05;
surf(C,D,temp,'edgecolor','none')
title('HARDI-ZOOMit Non-Interp','FontSize',14)
xlabel('parameter label','FontSize',14)
ylabel('parameter label','FontSize',14)
set(gca,'FontSize',12,...
        'ytick',1.5:(size(RHZni ,1)+0.5),...
        'yticklabel',table_names,...
        'xtick',1.5:(size(RHZni ,1)+0.5),...
        'xticklabel',table_names)
ytickangle(270)
colormap hot
colorbar
axis([ 1 size(RHZni ,1)+1 1 size(RHZni ,1)+1 min(RHZni (:)) max(RHZni (:)) ])
az = 90;
el = 90;
view(az, el);
