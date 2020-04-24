%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
clc
clear all
close all
%% dMRI parameter settings
save_path='/home/user/results';
workspace_file = fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);
%% Indexes and lsit of variables significant between-group differences for 60 dMRI acquisitions used in Labounek et al. (2020) Scientific Reports
var_indxs = [3 49 23 43 51 12 11 13 33 44 16 15 46 36 35 20 19 40 39 70 71 72 73 74 75];
table_names = {'FAwS' 'FAwK' 'f1wS' 'f1wSK' 'f1wK' 'MDwM' 'MDwm' 'MDwS' 'dwS' 'dwSK' ...
        'MDgM' 'MDgm' 'MDgSK' 'dgM' 'dgm' 'MDwgM' 'MDwgm' 'dwgM' 'dwgm' ...
        'FAwH' 'f1wH' 'MDwH' 'dwH' 'MDgH' 'dgH'};
table_pvals_rows = [9 63 28 57 65 20 45 23 33 61 21 46 60 31 51 22 49 32 52 75 76 77 78 81 82];
%% New protocol-specific group indexes
group1 = zeros(size(sbj));
group1(sbj==3) = 1;
group1(sbj==1) = 2;
group1(sbj==2) = 3;
group1(sbj==4) = 4;

group2 = zeros(size(sbj));
group2(sbj==3) = 5;
group2(sbj==1) = 6;
group2(sbj==2) = 7;
group2(sbj==4) = 8;
group2(sum(table_ZOOMit_NotInt,2)==0) = [];

group3 = zeros(size(sbj));
group3(sbj==3) = 9;
group3(sbj==1) = 10;
group3(sbj==2) = 11;
group3(sbj==4) = 12;

group = [group1; group2; group3]';
%% Figure definition
sizefig = [30 30 1300 870]; % for 4 rows 6 columns
% sizefig = [30 30 1300 1270];
space = [0.2 0.2];  %vertical, horizontal
rows = 4;
cols = 6;
pos = [0.05 0.05 0.9 0.9];

BigAx = newplot();
fig = ancestor(BigAx, 'figure');
set(fig, 'Position', sizefig);
hold_state = ishold(BigAx);
set(BigAx, 'Visible', 'off', 'color', 'none');

% pos = get(BigAx, 'Position');
set(BigAx, 'Position', pos);
width = pos(3)/cols;
height = pos(4)/rows; 

pos(1:2) = pos(1:2)+space.*[width height];

xlim = zeros([rows cols 2]);
ylim = zeros([rows cols 2]);

BigAxHV = get(BigAx, 'HandleVisibility');
BigAxParent = get(BigAx, 'Parent');
%% For loop drawing boxplots for each significant dMRI variable
ind = 1;
for i = 1:rows
    for j = 1:cols
            if ind<=25
                axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height ...
                    width*(1-space(1)) height*(1-space(2))];
                ax(i,j) = axes('Position', axPos, 'HandleVisibility', BigAxHV, 'parent', BigAxParent);
                set(ax(i,j),'visible', 'on');
                data1 = table_ZOOMit_Int(:,var_indxs(ind))';
                data2 = table_ZOOMit_NotInt(:,var_indxs(ind))';
                data2(sum(table_ZOOMit_NotInt,2)==0) = [];
                data3 = table_RESOLVE(:,var_indxs(ind))';
                data = [data1 data2 data3];
                if strcmp(table_names{1,ind},'MDwS')
                        ymax = 0.8;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'FAwS')
                        ymax = 0.22;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'FAwK')
                        ymax = 4.1;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'f1wS')
                        ymax = 0.22;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'f1wK')
                        ymax = 3.53;
                        triangle_pos = 0.91*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'MDwm')
                        ymax = 1.7;
                        triangle_pos = 0.91*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'MDwK')
                        ymax = 8.0;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'dwS')
                        ymax = 1.0;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'dwSK')
                        ymax=1.35*max(data);
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.90*ymax;
                elseif strcmp(table_names{1,ind},'MDgM') || strcmp(table_names{1,ind},'MDgm')
                        ymax=1.15*max(data);
                        triangle_pos = 0.90*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'MDgSK')
                        ymax=3.5;
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.90*ymax;
                elseif strcmp(table_names{1,ind},'dgM') || strcmp(table_names{1,ind},'dgm')
                         ymax=1.15*max(data);
                        triangle_pos = 0.90*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'f1wH') || strcmp(table_names{1,ind},'MDwH')
                        ymax=1.25*max(data);
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                elseif strcmp(table_names{1,ind},'f1wSK')
                        ymax=1.15*max(data);
                        triangle_pos = 0.75*ymax;
                        dot_pos = 0.90*ymax;
                else
                        ymax=1.15*max(data);
                        triangle_pos = 0.88*ymax;
                        dot_pos = 0.95*ymax;
                end
                ymin=min(data);
                if ymin <0
                        ymin = 1.1*ymin;
                else
                        ymin = 0.9*ymin;
                end
                will_pvals = table_pvals(table_pvals_rows(1,ind),[2 4 8 10 14 16 3 5 9 11 15 17]);
                draw_boxplot_grp(data,group,ymin,ymax,ind,will_pvals,WillThr,triangle_pos,dot_pos)
                title(table_names{1,ind})
                ind = ind+1;
            end
    end
end
% suptitle(supertitle)