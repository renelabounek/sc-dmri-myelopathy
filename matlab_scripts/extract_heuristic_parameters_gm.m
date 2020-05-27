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
% the Free Software Foundation, either version 3 of the License, or any later version.
%
% sc-dmri-myelopathy is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with sc-dmri-myelopathy.  If not, see <https://www.gnu.org/licenses/>.
%
clear all;
%% INITIALIZE AND DEFINE RESULT SOURCE AND SAVE FILE
% !!!! EXECUTION OF THE extract_descriptive_statistics.m AND extract_heuristic_parameters_gm.m SCRIPTS MUST
% PREVENT THIS SCRIPT CALL IN THE ORDER AS IT IS WRITTEN, OTHERWISE THE CURRENT SCRIPT WILL NOT WORK !!!!
save_path='/home/user/results';
workspace_file=fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);
% The original loaded workspace_file can have set different value,therefore it is defined here twice.
save_path='/home/user/results';
workspace_file=fullfile(save_path,'dmri_comparison_pvaltable.mat');
%% DOMAIN AND LIMITS INITIALIZATION
% Define dMRI metric domains for fill commands plotting confidece intervals.
FA_ax = [0, FAx_vec, FAx_vec(end:-1:1)];
f1_ax = [0, f1x_vec, f1x_vec(end:-1:1)];
MD_ax = [0, MDx_vec, MDx_vec(end:-1:1)];
d_ax = [0, dx_vec, dx_vec(end:-1:1)];

% Define positions of heuristic parameter limits
FA_heu_min = find(FAx_vec>=0.47, 1 );
FA_heu_max = find(FAx_vec<=0.67, 1, 'last' );
f1_heu_min = find(f1x_vec>=0.30, 1 );
f1_heu_max = find(f1x_vec<=0.55, 1, 'last' );
MD_heu_min = find(MDx_vec>=0.84, 1 );
MD_heu_max = find(MDx_vec<=1.26, 1, 'last' );
% MD_heu_min = find(MDx_vec>=0.72, 1 );
% MD_heu_max = find(MDx_vec<=0.93, 1, 'last' );
d_heu_min = find(dx_vec>=1.00, 1 );
d_heu_max = find(dx_vec<=1.48, 1, 'last' );
% These lemits did not demonstrate statistically significant differences now,
% but may become useful in the future research. So they are kept in the code.
FA_heu2_min = find(FAx_vec>=0.75, 1 );
FA_heu2_max = find(FAx_vec<=0.90, 1, 'last' );
f1_heu2_min = find(f1x_vec>=0.60, 1 );
f1_heu2_max = find(f1x_vec<=0.80, 1, 'last' );
%% GM HEURISTIC PARAMETER ESTIMATION AND DMRI METRIC SMOOTH DISTRIBUTION VISUALIZATION 
% Initialize figure 1.
h(1).fig = figure(1);
set(h(1).fig, 'Position', [50, 50, 1200, 950]);

% Deviation of FA GM heuristic parameters and their group-based
% distribution plots, i.e. Fig. 3 in Labounek et al. (2020) Scientific Reports
heu_sbj = [];
heu_prot = []; 
heu_FA = [];
heu_FA2 = [];
for prot = 1:size(FA_gm_kernel,3)
    pdfs = FA_gm_kernel(:,:,prot);
    pdfs_sum = sum(pdfs);
    pdfs_sbj = sbj;
    pdfs_sbj(pdfs_sum == 0) = [];
    pdfs(:,pdfs_sum == 0) = [];
    
    control_mean = mean(pdfs(:,pdfs_sbj==1),2);
    control_median = median(pdfs(:,pdfs_sbj==1),2);
    control_std = std(pdfs(:,pdfs_sbj==1),0,2);
    asympatlight_mean = mean(pdfs(:,pdfs_sbj==2),2);
    asympatlight_std = std(pdfs(:,pdfs_sbj==2),0,2);
    repro_mean = mean(pdfs(:,pdfs_sbj==3),2);
    repro_std = std(pdfs(:,pdfs_sbj==3),0,2);
    asympatstrong_mean = mean(pdfs(:,pdfs_sbj==4),2);
    asympatstrong_std = std(pdfs(:,pdfs_sbj==4),0,2);
    
    heu_FA_prot = sum(pdfs(FA_heu_min:FA_heu_max,:))';
    heu_FA2_prot = sum(pdfs(FA_heu2_min:FA_heu2_max,:))';
    heu_FA = [ heu_FA; heu_FA_prot];
    heu_FA2 = [ heu_FA2; heu_FA2_prot];
    heu_sbj = [heu_sbj; pdfs_sbj];
    heu_prot = [heu_prot; prot*ones(size(pdfs_sbj))];
    heu_FA_mean(1,prot) = mean(heu_FA_prot);
    heu_FA_median(1,prot) = median(heu_FA_prot);
    heu_FA2_mean(1,prot) = mean(heu_FA2_prot);
    heu_FA2_median(1,prot) = median(heu_FA2_prot);

    pom = pdfs(:,pdfs_sbj==1);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    control_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==3);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    repro_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==2);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatlight_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==4);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatstrong_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    subplot(4,3,prot)
    fill(FA_ax, repro_interval, 1,....
        'facecolor',0.9*[0.1 0.1 0.1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    hold on
%     plot(FAx_vec,control_median,'--','Color',0.8*[1 0 0], 'LineWidth', 2)
    fill(FA_ax, control_interval, 1,....
        'facecolor',0.9*[1 0 0], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(FA_ax, asympatlight_interval, 1,....
        'facecolor',0.9*[0 0.749 1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(FA_ax, asympatstrong_interval, 1,....
        'facecolor',0.9*[0.65 0.65 0.65], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    plot(FAx_vec,repro_mean,'Color',0.8*[0 0 0], 'LineWidth', 3)
    plot(FAx_vec,control_mean,'Color',0.8*[1 0 0], 'LineWidth', 3)
    plot(FAx_vec,asympatlight_mean,'Color',0.8*[0 0.749 1], 'LineWidth', 3)
    plot(FAx_vec,asympatstrong_mean,'Color',0.8*[0.65 0.65 0.65], 'LineWidth', 3)
    hold off
    if prot == 1
        ylabel('FA probabilities','FontSize',16)
        title('HARDI-ZOOMit Interp','FontSize',20)
    elseif prot == 2
        title('HARDI-ZOOMit Non-Interp','FontSize',20)
    elseif prot == 3
        title('DTI-RESOLVE Non-Interp','FontSize',20)
    end
    set(gca,'FontSize',14)
    axis([0 1 0 1.05*max([repro_interval; control_interval; asympatlight_interval])])
    pause(.1)
end

% Deviation of f1 GM heuristic parameters and their group-based
% distribution plots, i.e. Fig. 3 in Labounek et al. (2020) Scientific Reports
heu_f1 = [];
heu_f12 = [];
for prot = 1:size(f1_gm_kernel,3)
    pdfs = f1_gm_kernel(:,:,prot);
    pdfs_sum = sum(pdfs);
    pdfs_sbj = sbj;
    pdfs_sbj(pdfs_sum == 0) = [];
    pdfs(:,pdfs_sum == 0) = [];
    
    control_mean = mean(pdfs(:,pdfs_sbj==1),2);
    control_median = median(pdfs(:,pdfs_sbj==1),2);
    control_std = std(pdfs(:,pdfs_sbj==1),0,2);
    asympatlight_mean = mean(pdfs(:,pdfs_sbj==2),2);
    asympatlight_std = std(pdfs(:,pdfs_sbj==2),0,2);
    repro_mean = mean(pdfs(:,pdfs_sbj==3),2);
    repro_std = std(pdfs(:,pdfs_sbj==3),0,2);
    asympatstrong_mean = mean(pdfs(:,pdfs_sbj==4),2);
    asympatstrong_std = std(pdfs(:,pdfs_sbj==4),0,2);
    
    heu_f1_prot = sum(pdfs(f1_heu_min:f1_heu_max,:))';
    heu_f12_prot = sum(pdfs(f1_heu2_min:f1_heu2_max,:))';
    heu_f1 = [ heu_f1; heu_f1_prot ];
    heu_f12 = [ heu_f12; heu_f12_prot ];
    heu_f1_mean(1,prot) = mean(heu_f1_prot);
    heu_f1_median(1,prot) = median(heu_f1_prot);
    heu_f12_mean(1,prot) = mean(heu_f12_prot);
    heu_f12_median(1,prot) = median(heu_f12_prot);
    
    pom = pdfs(:,pdfs_sbj==1);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    control_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==3);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    repro_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==2);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatlight_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==4);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatstrong_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    subplot(4,3,3+prot)
    fill(f1_ax, repro_interval, 1,....
        'facecolor',0.9*[0.1 0.1 0.1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    hold on
%     plot(f1x_vec,control_median,'--','Color',0.8*[1 0 0], 'LineWidth', 2)
    fill(f1_ax, control_interval, 1,....
        'facecolor',0.9*[1 0 0], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(f1_ax, asympatlight_interval, 1,....
        'facecolor',0.9*[0 0.749 1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(f1_ax, asympatstrong_interval, 1,....
        'facecolor',0.9*[0.65 0.65 0.65], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    plot(f1x_vec,repro_mean,'Color',0.8*[0 0 0], 'LineWidth', 3)
    plot(f1x_vec,control_mean,'Color',0.8*[1 0 0], 'LineWidth', 3)
    plot(f1x_vec,asympatlight_mean,'Color',0.8*[0 0.749 1], 'LineWidth', 3)
    plot(f1x_vec,asympatstrong_mean,'Color',0.8*[0.65 0.65 0.65], 'LineWidth', 3)
    hold off
    if prot == 1
        ylabel('f_1 probabilities','FontSize',16)
    end
    set(gca,'FontSize',14)
    axis([0 1 0 1.05*max([repro_interval; control_interval; asympatlight_interval])])
    pause(.1)
end

% Deviation of d GM heuristic parameters and their group-based
% distribution plots, i.e. Fig. 3 in Labounek et al. (2020) Scientific Reports
heu_d = [];
for prot = 1:size(d_gm_kernel,3)
    pdfs = d_gm_kernel(:,:,prot);
    pdfs_sum = sum(pdfs);
    pdfs_sbj = sbj;
    pdfs_sbj(pdfs_sum == 0) = [];
    pdfs(:,pdfs_sum == 0) = [];
    indxs = find(pdfs_sum ~= 0);
    
    control_mean = mean(pdfs(:,pdfs_sbj==1),2);
    control_median = median(pdfs(:,pdfs_sbj==1),2);
    control_std = std(pdfs(:,pdfs_sbj==1),0,2);
    asympatlight_mean = mean(pdfs(:,pdfs_sbj==2),2);
    asympatlight_std = std(pdfs(:,pdfs_sbj==2),0,2);
    repro_mean = mean(pdfs(:,pdfs_sbj==3),2);
    repro_std = std(pdfs(:,pdfs_sbj==3),0,2);
    asympatstrong_mean = mean(pdfs(:,pdfs_sbj==4),2);
    asympatstrong_std = std(pdfs(:,pdfs_sbj==4),0,2);
    
    heu_d_prot = sum(pdfs(d_heu_min:d_heu_max,:))';
    heu_d = [ heu_d; heu_d_prot ];
    heu_d_mean(1,prot) = mean(heu_d_prot);
    heu_d_median(1,prot) = median(heu_d_prot);
    
    table_vec = zeros(size(FA_wm_mean,1),1);
    table_vec(indxs,1) = heu_d_prot;
    
    if prot == 1
        table_ZOOMit_Int(:,75) = table_vec;
    elseif prot == 2
        table_ZOOMit_NotInt(:,75) = table_vec;
    elseif prot == 3
        table_RESOLVE(:,75) = table_vec;
    end
    
    pom = pdfs(:,pdfs_sbj==1);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    control_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==3);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    repro_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==2);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatlight_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==4);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatstrong_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    subplot(4,3,9+prot)
    fill(d_ax, repro_interval, 1,....
        'facecolor',0.9*[0.1 0.1 0.1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    hold on
%     plot(dx_vec,control_median,'--','Color',0.8*[1 0 0], 'LineWidth', 2)
    fill(d_ax, control_interval, 1,....
        'facecolor',0.9*[1 0 0], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(d_ax, asympatlight_interval, 1,....
        'facecolor',0.9*[0 0.749 1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(d_ax, asympatstrong_interval, 1,....
        'facecolor',0.9*[0.65 0.65 0.65], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    H1 = plot(dx_vec,repro_mean,'Color',0.8*[0 0 0], 'LineWidth', 3);
    H2 = plot(dx_vec,control_mean,'Color',0.8*[1 0 0], 'LineWidth', 3);
    H3 = plot(dx_vec,asympatlight_mean,'Color',0.8*[0 0.749 1], 'LineWidth',3);
    H4 =  plot(dx_vec,asympatstrong_mean,'Color',0.8*[0.65 0.65 0.65], 'LineWidth', 3);
    hold off
    if prot == 1
        ylabel('d probabilities','FontSize',16)       
    end
    set(gca,'FontSize',14)
    xlabel('*10^{-9}m^2/s','FontSize',13)
    axis([0 4 0 1.05*max([repro_interval; control_interval; asympatlight_interval])])
    pause(.1)
end

% Deviation of MD GM heuristic parameters and their group-based
% distribution plots, i.e. Fig. 3 in Labounek et al. (2020) Scientific Reports
heu_MD = [];
for prot = 1:size(MD_gm_kernel,3)
    pdfs = MD_gm_kernel(:,:,prot);
    pdfs_sum = sum(pdfs);
    pdfs_sbj = sbj;
    pdfs_sbj(pdfs_sum == 0) = [];
    pdfs(:,pdfs_sum == 0) = [];
    indxs = find(pdfs_sum ~= 0);
    
    control_mean = mean(pdfs(:,pdfs_sbj==1),2);
    control_median = median(pdfs(:,pdfs_sbj==1),2);
    control_std = std(pdfs(:,pdfs_sbj==1),0,2);
    asympatlight_mean = mean(pdfs(:,pdfs_sbj==2),2);
    asympatlight_std = std(pdfs(:,pdfs_sbj==2),0,2);
    repro_mean = mean(pdfs(:,pdfs_sbj==3),2);
    repro_std = std(pdfs(:,pdfs_sbj==3),0,2);
    asympatstrong_mean = mean(pdfs(:,pdfs_sbj==4),2);
    asympatstrong_std = std(pdfs(:,pdfs_sbj==4),0,2);
    
    heu_MD_prot = sum(pdfs(MD_heu_min:MD_heu_max,:))';
    heu_MD = [ heu_MD; heu_MD_prot ];
    heu_MD_mean(1,prot) = mean(heu_MD_prot);
    heu_MD_median(1,prot) = median(heu_MD_prot);
    
    table_vec = zeros(size(FA_wm_mean,1),1);
    table_vec(indxs,1) = heu_MD_prot;
    
    if prot == 1
        table_ZOOMit_Int(:,74) = table_vec;
    elseif prot == 2
        table_ZOOMit_NotInt(:,74) = table_vec;
    elseif prot == 3
        table_RESOLVE(:,74) = table_vec;
    end

    pom = pdfs(:,pdfs_sbj==1);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    control_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==3);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    repro_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==2);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatlight_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    pom = pdfs(:,pdfs_sbj==4);
    for t = 1:size(pom,1)
        y = sort(pom(t,:));
        % compute 25th percentile (first quartile)
        Q(t,1) = median(y(find(y<median(y))));
        % compute 50th percentile (second quartile)
        Q(t,2) = median(y);
        % compute 75th percentile (third quartile)
        Q(t,3) = median(y(find(y>median(y))));
    end
    asympatstrong_interval = [Q(1,3); Q(:,1); Q(end:-1:1,3)];
    
    subplot(4,3,6+prot)
    fill(MD_ax, repro_interval, 1,....
        'facecolor',0.9*[0.1 0.1 0.1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    hold on
%     plot(MDx_vec,control_median,'--','Color',0.8*[1 0 0], 'LineWidth', 2)
    fill(MD_ax, control_interval, 1,....
        'facecolor',0.9*[1 0 0], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(MD_ax, asympatlight_interval, 1,....
        'facecolor',0.9*[0 0.749 1], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    fill(MD_ax, asympatstrong_interval, 1,....
        'facecolor',0.9*[0.65 0.65 0.65], ...
        'edgecolor','none', ...
        'facealpha', 0.25);
    H1 = plot(MDx_vec,repro_mean,'Color',0.8*[0 0 0], 'LineWidth', 3);
    H2 = plot(MDx_vec,control_mean,'Color',0.8*[1 0 0], 'LineWidth', 3);
    H3 = plot(MDx_vec,asympatlight_mean,'Color',0.8*[0 0.749 1], 'LineWidth', 3);
    H4 = plot(MDx_vec,asympatstrong_mean,'Color',0.8*[0.65 0.65 0.65], 'LineWidth', 3);
    plot([3 4],[0 0],'Color',0.8*[0.65 0.65 0.65], 'LineWidth', 3)
    hold off
    if prot == 1
        ylabel('MD probabilities','FontSize',16)
    elseif prot == 3
%         xlabel('RESOLVE','FontSize',16)
        legend([H1 H2 H3 H4],{'Reproducibility','Control','MC Asym. Patients','SC Asym. Patients'})
        set(h(1).fig.Children(1,1),'units','pixel')
        set(h(1).fig.Children(1,1),'FontSize',12)
        set(h(1).fig.Children(1,1),'position',[1045 200 110 70])    
    end
    set(gca,'FontSize',14)
    xlabel('*10^{-9}m^2/s','FontSize',13)
    axis([0 4 0 1.05*max([repro_interval; control_interval; asympatlight_interval])])
    pause(.1)
end
pause(3)
print(fullfile(save_path,'dmri_distributions_gm'), '-dpng', '-r300')
close(h(1).fig)
pause(1)

%% GROUP DESCRIPTIVE STATISTICS FOR ALL DERIVED DMRI METRICS STORED IN PROTOCOL-SPECIFIC TABLES
for category = 1:3
    par_mat = table_ZOOMit_Int(sbj==category,:);
    par_mat(sum(par_mat,2)==0,:) = [];
    table_ZOOMit_Int_mean(category,1:75) = mean(par_mat);
    table_ZOOMit_Int_median(category,1:75) = median(par_mat);
    table_ZOOMit_Int_std(category,1:75) = std(par_mat);

    par_mat = table_ZOOMit_NotInt(sbj==category,:);
    par_mat(sum(par_mat,2)==0,:) = [];
    table_ZOOMit_NotInt_mean(category,1:75) = mean(par_mat);
    table_ZOOMit_NotInt_median(category,1:75) = median(par_mat);
    table_ZOOMit_NotInt_std(category,1:75) = std(par_mat);

    par_mat = table_RESOLVE(sbj==category,:);
    par_mat(sum(par_mat,2)==0,:) = [];
    table_RESOLVE_mean(category,1:75) = mean(par_mat);
    table_RESOLVE_median(category,1:75) = median(par_mat);
    table_RESOLVE_std(category,1:75) = std(par_mat);
end

%% STORING WM HEURISTIC RESULTS INTO GM SPECIFIC HEU VARIABLES
heu_MDGM = heu_MD;
heu_dGM = heu_d;
heu_FAGM = heu_FA;
heu_f1GM = heu_f1;

%% VISUALIZATION OF FA HEURISTIC PARAMETER DISTRIBUTIONS FROM GM
pCR_ZOOMint = ranksum(heu_FA(heu_sbj==1 & heu_prot==1),heu_FA(heu_sbj==3 & heu_prot==1));
pCP_ZOOMint = ranksum(heu_FA(heu_sbj==1 & heu_prot==1),heu_FA(heu_sbj==2 & heu_prot==1));
pRP_ZOOMint = ranksum(heu_FA(heu_sbj==3 & heu_prot==1),heu_FA(heu_sbj==2 & heu_prot==1));
pSC_ZOOMint = ranksum(heu_FA(heu_sbj==4 & heu_prot==1),heu_FA(heu_sbj==1 & heu_prot==1));
pSR_ZOOMint = ranksum(heu_FA(heu_sbj==3 & heu_prot==1),heu_FA(heu_sbj==4 & heu_prot==1));
pSP_ZOOMint = ranksum(heu_FA(heu_sbj==4 & heu_prot==1),heu_FA(heu_sbj==2 & heu_prot==1));
pCR_ZOOMnotint = ranksum(heu_FA(heu_sbj==1 & heu_prot==2),heu_FA(heu_sbj==3 & heu_prot==2));
pCP_ZOOMnotint = ranksum(heu_FA(heu_sbj==1 & heu_prot==2),heu_FA(heu_sbj==2 & heu_prot==2));
pRP_ZOOMnotint = ranksum(heu_FA(heu_sbj==3 & heu_prot==2),heu_FA(heu_sbj==2 & heu_prot==2));
pSC_ZOOMnotint = ranksum(heu_FA(heu_sbj==4 & heu_prot==2),heu_FA(heu_sbj==1 & heu_prot==2));
pSR_ZOOMnotint = ranksum(heu_FA(heu_sbj==3 & heu_prot==2),heu_FA(heu_sbj==4 & heu_prot==2));
pSP_ZOOMnotint = ranksum(heu_FA(heu_sbj==4 & heu_prot==2),heu_FA(heu_sbj==2 & heu_prot==2));
pCR_RESOLVE = ranksum(heu_FA(heu_sbj==1 & heu_prot==3),heu_FA(heu_sbj==3 & heu_prot==3));
pCP_RESOLVE = ranksum(heu_FA(heu_sbj==1 & heu_prot==3),heu_FA(heu_sbj==2 & heu_prot==3));
pRP_RESOLVE = ranksum(heu_FA(heu_sbj==3 & heu_prot==3),heu_FA(heu_sbj==2 & heu_prot==3));
pSC_RESOLVE = ranksum(heu_FA(heu_sbj==4 & heu_prot==3),heu_FA(heu_sbj==1 & heu_prot==3));
pSR_RESOLVE = ranksum(heu_FA(heu_sbj==3 & heu_prot==3),heu_FA(heu_sbj==4 & heu_prot==3));
pSP_RESOLVE = ranksum(heu_FA(heu_sbj==4 & heu_prot==3),heu_FA(heu_sbj==2 & heu_prot==3));
table_pvals(79,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
h(2).fig = figure(2);
set(h(2).fig, 'Position', [50, 50, 750, 550]);
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_FA_mean',1,2)','c-','LineWidth',5)
hold on
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_FA_median',1,2)','m-','LineWidth',5)
scatter(heu_prot(heu_sbj==3)/2, heu_FA(heu_sbj==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.7);
scatter(heu_prot(heu_sbj==1)/2, heu_FA(heu_sbj==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.7);
scatter(heu_prot(heu_sbj==4)/2, heu_FA(heu_sbj==4),850, '.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
scatter(heu_prot(heu_sbj==2)/2, heu_FA(heu_sbj==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.7);
if pCR_ZOOMint < WillThr
    plot((1-0.31)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((1-0.25)/2,0.785,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMint < WillThr
    plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMint < WillThr
    plot((1-0.31)/2,0.735,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.735,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMint < WillThr
    plot((1-0.31)/2,0.710,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.710,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMint < WillThr
    plot((1-0.31)/2,0.685,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.685,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMint < WillThr
    plot((1-0.31)/2,0.660,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.660,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((2-0.25)/2,0.785,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.735,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.735,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.710,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.710,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.685,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.685,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.660,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.660,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCR_RESOLVE < WillThr
    plot((3-0.37)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((3-0.31)/2,0.785,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pCP_RESOLVE < WillThr
    plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.760,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pRP_RESOLVE < WillThr
    plot((3-0.37)/2,0.735,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.735,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSC_RESOLVE < WillThr
    plot((3-0.37)/2,0.710,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.710,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSR_RESOLVE < WillThr
    plot((3-0.37)/2,0.685,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.685,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSP_RESOLVE < WillThr
    plot((3-0.37)/2,0.660,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.660,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
hold off
title({'FA heuristic parameter probability distribution from GM';'C3-C6 segments over subjects'})
set(gca,'XTick',(1:1:3)/2,...
     'XTickLabel',{'ZOOMit interp'
                   'ZOOMit non-interp'
                   'RESOLVE'
                   },...
     'TickLength',[0 0],'LineWidth',2,...
     'FontSize',14)
xlabel('Different protocols','FontSize',18)
ylabel('Probability of FA heuristic interval','FontSize',18)
axis([0.5/2 3.5/2 0.2 0.80])
% axis([0.25 1.75 0.2 0.8])
grid on
pause(1)
print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_heu_GM_over_subjects']), '-dpng', '-r300')

%% VISUALIZATION OF f1 HEURISTIC PARAMETER DISTRIBUTIONS FROM GM
pCR_ZOOMint = ranksum(heu_f1(heu_sbj==1 & heu_prot==1),heu_f1(heu_sbj==3 & heu_prot==1));
pCP_ZOOMint = ranksum(heu_f1(heu_sbj==1 & heu_prot==1),heu_f1(heu_sbj==2 & heu_prot==1));
pRP_ZOOMint = ranksum(heu_f1(heu_sbj==3 & heu_prot==1),heu_f1(heu_sbj==2 & heu_prot==1));
pSC_ZOOMint = ranksum(heu_f1(heu_sbj==4 & heu_prot==1),heu_f1(heu_sbj==1 & heu_prot==1));
pSR_ZOOMint = ranksum(heu_f1(heu_sbj==3 & heu_prot==1),heu_f1(heu_sbj==4 & heu_prot==1));
pSP_ZOOMint = ranksum(heu_f1(heu_sbj==4 & heu_prot==1),heu_f1(heu_sbj==2 & heu_prot==1));
pCR_ZOOMnotint = ranksum(heu_f1(heu_sbj==1 & heu_prot==2),heu_f1(heu_sbj==3 & heu_prot==2));
pCP_ZOOMnotint = ranksum(heu_f1(heu_sbj==1 & heu_prot==2),heu_f1(heu_sbj==2 & heu_prot==2));
pRP_ZOOMnotint = ranksum(heu_f1(heu_sbj==3 & heu_prot==2),heu_f1(heu_sbj==2 & heu_prot==2));
pSC_ZOOMnotint = ranksum(heu_f1(heu_sbj==4 & heu_prot==2),heu_f1(heu_sbj==1 & heu_prot==2));
pSR_ZOOMnotint = ranksum(heu_f1(heu_sbj==3 & heu_prot==2),heu_f1(heu_sbj==4 & heu_prot==2));
pSP_ZOOMnotint = ranksum(heu_f1(heu_sbj==4 & heu_prot==2),heu_f1(heu_sbj==2 & heu_prot==2));
pCR_RESOLVE = ranksum(heu_f1(heu_sbj==1 & heu_prot==3),heu_f1(heu_sbj==3 & heu_prot==3));
pCP_RESOLVE = ranksum(heu_f1(heu_sbj==1 & heu_prot==3),heu_f1(heu_sbj==2 & heu_prot==3));
pRP_RESOLVE = ranksum(heu_f1(heu_sbj==3 & heu_prot==3),heu_f1(heu_sbj==2 & heu_prot==3));
pSC_RESOLVE = ranksum(heu_f1(heu_sbj==4 & heu_prot==3),heu_f1(heu_sbj==1 & heu_prot==3));
pSR_RESOLVE = ranksum(heu_f1(heu_sbj==3 & heu_prot==3),heu_f1(heu_sbj==4 & heu_prot==3));
pSP_RESOLVE = ranksum(heu_f1(heu_sbj==4 & heu_prot==3),heu_f1(heu_sbj==2 & heu_prot==3));
table_pvals(80,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
h(3).fig = figure(3);
set(h(3).fig, 'Position', [50, 50, 750, 550]);
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_f1_mean',1,2)','c-','LineWidth',5)
hold on
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_f1_median',1,2)','m-','LineWidth',5)
scatter(heu_prot(heu_sbj==3)/2, heu_f1(heu_sbj==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.8);
scatter(heu_prot(heu_sbj==1)/2, heu_f1(heu_sbj==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.8);
scatter(heu_prot(heu_sbj==4)/2, heu_f1(heu_sbj==4),850, '.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeColor',[0.75 0.75 0.75],'MarkerEdgeAlpha',0.8);
scatter(heu_prot(heu_sbj==2)/2, heu_f1(heu_sbj==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.8);
if pCR_ZOOMint < WillThr
    plot((1-0.31)/2,0.885,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((1-0.25)/2,0.885,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMint < WillThr
    plot((1-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.860,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMint < WillThr
    plot((1-0.31)/2,0.835,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.835,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMint < WillThr
    plot((1-0.31)/2,0.810,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.810,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMint < WillThr
    plot((1-0.31)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.785,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMint < WillThr
    plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.760,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.885,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((2-0.25)/2,0.885,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.860,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.835,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.835,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.810,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.810,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.785,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.760,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCR_RESOLVE < WillThr
    plot((3-0.37)/2,0.885,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((3-0.31)/2,0.885,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pCP_RESOLVE < WillThr
    plot((3-0.37)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.860,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pRP_RESOLVE < WillThr
    plot((3-0.37)/2,0.835,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.835,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSC_RESOLVE < WillThr
    plot((3-0.37)/2,0.810,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.810,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSR_RESOLVE < WillThr
    plot((3-0.37)/2,0.785,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.785,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSP_RESOLVE < WillThr
    plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.760,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
hold off
title({'f_1 heuristic parameter probability distribution from GM';'C3-C6 segments over subjects'})
set(gca,'XTick',(1:1:3)/2,...
     'XTickLabel',{'ZOOMit interp'
                   'ZOOMit non-interp'
                   'RESOLVE'
                   },...
     'TickLength',[0 0],'LineWidth',2,...
     'FontSize',14)
xlabel('Different protocols','FontSize',18)
ylabel('Probability of f_1 heuristic interval','FontSize',18)
axis([0.5/2 3.5/2 0.2 0.90])
% axis([0.25 1.75 0.2 0.8])
grid on
pause(1)
print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_heu_GM_over_subjects']), '-dpng', '-r300')

%% VISUALIZATION OF MD HEURISTIC PARAMETER DISTRIBUTIONS FROM GM
pCR_ZOOMint = ranksum(heu_MD(heu_sbj==1 & heu_prot==1),heu_MD(heu_sbj==3 & heu_prot==1));
pCP_ZOOMint = ranksum(heu_MD(heu_sbj==1 & heu_prot==1),heu_MD(heu_sbj==2 & heu_prot==1));
pRP_ZOOMint = ranksum(heu_MD(heu_sbj==3 & heu_prot==1),heu_MD(heu_sbj==2 & heu_prot==1));
pSC_ZOOMint = ranksum(heu_MD(heu_sbj==4 & heu_prot==1),heu_MD(heu_sbj==1 & heu_prot==1));
pSR_ZOOMint = ranksum(heu_MD(heu_sbj==3 & heu_prot==1),heu_MD(heu_sbj==4 & heu_prot==1));
pSP_ZOOMint = ranksum(heu_MD(heu_sbj==4 & heu_prot==1),heu_MD(heu_sbj==2 & heu_prot==1));
pCR_ZOOMnotint = ranksum(heu_MD(heu_sbj==1 & heu_prot==2),heu_MD(heu_sbj==3 & heu_prot==2));
pCP_ZOOMnotint = ranksum(heu_MD(heu_sbj==1 & heu_prot==2),heu_MD(heu_sbj==2 & heu_prot==2));
pRP_ZOOMnotint = ranksum(heu_MD(heu_sbj==3 & heu_prot==2),heu_MD(heu_sbj==2 & heu_prot==2));
pSC_ZOOMnotint = ranksum(heu_MD(heu_sbj==4 & heu_prot==2),heu_MD(heu_sbj==1 & heu_prot==2));
pSR_ZOOMnotint = ranksum(heu_MD(heu_sbj==3 & heu_prot==2),heu_MD(heu_sbj==4 & heu_prot==2));
pSP_ZOOMnotint = ranksum(heu_MD(heu_sbj==4 & heu_prot==2),heu_MD(heu_sbj==2 & heu_prot==2));
pCR_RESOLVE = ranksum(heu_MD(heu_sbj==1 & heu_prot==3),heu_MD(heu_sbj==3 & heu_prot==3));
pCP_RESOLVE = ranksum(heu_MD(heu_sbj==1 & heu_prot==3),heu_MD(heu_sbj==2 & heu_prot==3));
pRP_RESOLVE = ranksum(heu_MD(heu_sbj==3 & heu_prot==3),heu_MD(heu_sbj==2 & heu_prot==3));
pSC_RESOLVE = ranksum(heu_MD(heu_sbj==4 & heu_prot==3),heu_MD(heu_sbj==1 & heu_prot==3));
pSR_RESOLVE = ranksum(heu_MD(heu_sbj==3 & heu_prot==3),heu_MD(heu_sbj==4 & heu_prot==3));
pSP_RESOLVE = ranksum(heu_MD(heu_sbj==4 & heu_prot==3),heu_MD(heu_sbj==2 & heu_prot==3));
table_pvals(81,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
h(4).fig = figure(4);
set(h(4).fig, 'Position', [50, 50, 750, 550]);
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_MD_mean',1,2)','c-','LineWidth',5)
hold on
plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_MD_median',1,2)','m-','LineWidth',5)
scatter(heu_prot(heu_sbj==3)/2, heu_MD(heu_sbj==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
scatter(heu_prot(heu_sbj==1)/2, heu_MD(heu_sbj==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
scatter(heu_prot(heu_sbj==4)/2, heu_MD(heu_sbj==4),850, '.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeColor',[0.85 0.85 0.85],'MarkerEdgeAlpha',0.9);
scatter(heu_prot(heu_sbj==2)/2, heu_MD(heu_sbj==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
if pCR_ZOOMint < WillThr
    plot((1-0.31)/2,1.070,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((1-0.25)/2,1.070,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMint < WillThr
    plot((1-0.31)/2,1.040,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,1.040,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMint < WillThr
    plot((1-0.31)/2,1.010,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,1.010,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMint < WillThr
    plot((1-0.31)/2,0.980,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.980,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMint < WillThr
    plot((1-0.31)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.950,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMint < WillThr
    plot((1-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.860,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCR_ZOOMnotint < WillThr
    plot((2-0.31)/2,1.070,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((2-0.25)/2,1.070,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMnotint < WillThr
    plot((2-0.31)/2,1.040,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,1.040,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMnotint < WillThr
    plot((2-0.31)/2,1.010,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,1.010,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.980,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.980,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.950,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.860,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCR_RESOLVE < WillThr
    plot((3-0.37)/2,1.070,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((3-0.31)/2,1.070,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pCP_RESOLVE < WillThr
    plot((3-0.37)/2,1.040,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,1.040,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pRP_RESOLVE < WillThr
    plot((3-0.37)/2,1.010,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,1.010,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSC_RESOLVE < WillThr
    plot((3-0.37)/2,0.980,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.980,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSR_RESOLVE < WillThr
    plot((3-0.37)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.950,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSP_RESOLVE < WillThr
    plot((3-0.37)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.860,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
hold off
title({'MD heuristic parameter probability distribution from GM';'C3-C6 segments over subjects'})
set(gca,'XTick',(1:1:3)/2,...
     'XTickLabel',{'ZOOMit interp'
                   'ZOOMit non-interp'
                   'RESOLVE'
                   },...
     'TickLength',[0 0],'LineWidth',2,...
     'FontSize',14)
xlabel('Different protocols','FontSize',18)
ylabel('Probability of MD heuristic interval','FontSize',18)
axis([0.5/2 3.5/2 0.2 1.09])
% axis([0.25 1.75 0.2 0.9])
grid on
pause(1)
print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_heu_GM_over_subjects']), '-dpng', '-r300')

%% VISUALIZATION OF d HEURISTIC PARAMETER DISTRIBUTIONS FROM GM
pCR_ZOOMint = ranksum(heu_d(heu_sbj==1 & heu_prot==1),heu_d(heu_sbj==3 & heu_prot==1));
pCP_ZOOMint = ranksum(heu_d(heu_sbj==1 & heu_prot==1),heu_d(heu_sbj==2 & heu_prot==1));
pRP_ZOOMint = ranksum(heu_d(heu_sbj==3 & heu_prot==1),heu_d(heu_sbj==2 & heu_prot==1));
pSC_ZOOMint = ranksum(heu_d(heu_sbj==4 & heu_prot==1),heu_d(heu_sbj==1 & heu_prot==1));
pSR_ZOOMint = ranksum(heu_d(heu_sbj==3 & heu_prot==1),heu_d(heu_sbj==4 & heu_prot==1));
pSP_ZOOMint = ranksum(heu_d(heu_sbj==4 & heu_prot==1),heu_d(heu_sbj==2 & heu_prot==1));
pCR_ZOOMnotint = ranksum(heu_d(heu_sbj==1 & heu_prot==2),heu_d(heu_sbj==3 & heu_prot==2));
pCP_ZOOMnotint = ranksum(heu_d(heu_sbj==1 & heu_prot==2),heu_d(heu_sbj==2 & heu_prot==2));
pRP_ZOOMnotint = ranksum(heu_d(heu_sbj==3 & heu_prot==2),heu_d(heu_sbj==2 & heu_prot==2));
pSC_ZOOMnotint = ranksum(heu_d(heu_sbj==4 & heu_prot==2),heu_d(heu_sbj==1 & heu_prot==2));
pSR_ZOOMnotint = ranksum(heu_d(heu_sbj==3 & heu_prot==2),heu_d(heu_sbj==4 & heu_prot==2));
pSP_ZOOMnotint = ranksum(heu_d(heu_sbj==4 & heu_prot==2),heu_d(heu_sbj==2 & heu_prot==2));
pCR_RESOLVE = ranksum(heu_d(heu_sbj==1 & heu_prot==3),heu_d(heu_sbj==3 & heu_prot==3));
pCP_RESOLVE = ranksum(heu_d(heu_sbj==1 & heu_prot==3),heu_d(heu_sbj==2 & heu_prot==3));
pRP_RESOLVE = ranksum(heu_d(heu_sbj==3 & heu_prot==3),heu_d(heu_sbj==2 & heu_prot==3));
pSC_RESOLVE = ranksum(heu_d(heu_sbj==4 & heu_prot==3),heu_d(heu_sbj==1 & heu_prot==3));
pSR_RESOLVE = ranksum(heu_d(heu_sbj==3 & heu_prot==3),heu_d(heu_sbj==4 & heu_prot==3));
pSP_RESOLVE = ranksum(heu_d(heu_sbj==4 & heu_prot==3),heu_d(heu_sbj==2 & heu_prot==3));
table_pvals(82,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
h(5).fig = figure(5);
set(h(5).fig, 'Position', [50, 50, 750, 550]);
HH1 = plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_d_mean',1,2)','c-','LineWidth',5);
hold on
HH2 = plot([unique(heu_prot)/2-0.15 unique(heu_prot)/2+0.15]', repmat(heu_d_median',1,2)','m-','LineWidth',5);
HH3 = scatter(heu_prot(heu_sbj==3)/2, heu_d(heu_sbj==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
HH4 = scatter(heu_prot(heu_sbj==1)/2, heu_d(heu_sbj==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
HH5 = scatter(heu_prot(heu_sbj==4)/2, heu_d(heu_sbj==4),850, '.', 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeColor',[0.85 0.85 0.85],'MarkerEdgeAlpha',0.9);
HH6 = scatter(heu_prot(heu_sbj==2)/2, heu_d(heu_sbj==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.12,'MarkerEdgeAlpha',0.9);
% legend([HH1(1,1) HH2(1,1) HH3 HH4 HH6 HH5],{'Whole dataset mean','Whole dataset median','Reproducibility',...
%     'Controls','Low comp. asym. patients','Strong comp. asym. patients'})
if pCR_ZOOMint < WillThr
    plot((1-0.31)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((1-0.25)/2,0.950,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMint < WillThr
    plot((1-0.31)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.900,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMint < WillThr
    plot((1-0.31)/2,0.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.850,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMint < WillThr
    plot((1-0.31)/2,0.910,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.910,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMint < WillThr
    plot((1-0.31)/2,0.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.800,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMint < WillThr
    plot((1-0.31)/2,0.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((1-0.25)/2,0.750,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
end
if pCR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((2-0.25)/2,0.950,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.900,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pRP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.850,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSC_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.910,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.910,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSR_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.800,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pSP_ZOOMnotint < WillThr
    plot((2-0.31)/2,0.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((2-0.25)/2,0.750,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
end
if pCR_RESOLVE < WillThr
    plot((3-0.37)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
    text((3-0.31)/2,0.950,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pCP_RESOLVE < WillThr
    plot((3-0.37)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.900,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pRP_RESOLVE < WillThr
    plot((3-0.37)/2,0.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.850,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSC_RESOLVE < WillThr
    plot((3-0.37)/2,0.980,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.980,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSR_RESOLVE < WillThr
    plot((3-0.37)/2,0.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.950,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
end
if pSP_RESOLVE < WillThr
    plot((3-0.37)/2,0.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
    text((3-0.31)/2,0.750,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
end
hold off
title({'d heuristic parameter probability distribution from GM';'C3-C6 segments over subjects'})
set(gca,'XTick',(1:1:3)/2,...
     'XTickLabel',{'ZOOMit interp'
                   'ZOOMit non-interp'
                   'RESOLVE'
                   },...
     'TickLength',[0 0],'LineWidth',2,...
     'FontSize',14)
xlabel('Different protocols','FontSize',18)
ylabel('Probability of d heuristic interval','FontSize',18)
axis([0.5/2 3.5/2 0.0 1.00])
% axis([0.25 1.75 0.2 0.9])
grid on
pause(1)
print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_heu_GM_over_subjects']), '-dpng', '-r300')

%% SAVE NEW GM HEURISTIC RESULTS INTO THE ORIGINAL workspace_file
% Now the file will contain all derived dmri-metrics and will be ready and stored for further analyses.
close all
pause(2)
save(workspace_file)
