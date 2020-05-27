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
%% INPUT PARAMETER SETTINGS
data_folder='/home/user/data'; %inside data folder there are subject specific subfoldres S0001, S0002,..., S00067
% as all uploded at URL: https://hdl.handle.net/20.500.12618/0000-5c13d342-4798-41d9-8d2a-bf750ab79fdb

% Your results will be stored in this folder:
save_path='/home/user/results';

% read_data equal to 1, if you want to read processed data and extract metrics from them
% read_data equal to 0, if you already have the extracted metrics and your wish is to draw graphs, i.e. you already have the workspace_file = fullfile(save_path,'dmri_comparison.mat')
read_data = 1;

% plot_results equal to 1, if you wish to plot graphs with dMRI metric distributions
% i.e. FigS 1,2,3,4,8 in Labounek et al. (2020) Scientific Reports
plot_results = 1;

% Set threshold for Wilcoxon rank-sum tests
WillThr = 0.05/6;

% Labounek et al. (2020) Scientific Reports subject list.
% The used subject list and their group categories (classified by radiologists and neurologists) were as follows:
% 1 - age-comparable helathy control, 2 - mild ADCCC patient, 3 - reproducibility group, 4 - severe ADCCC patient
if read_data == 1
    subject = {...
        'S0001'	2;...
        'S0002'	2;...
        'S0003'	1;...
        'S0004'	1;...
        'S0005'	3;...
        'S0006'	3;...
        'S0007'	4;...
        'S0008'	4;...
        'S0009'	3;...
        'S0010'	3;...
        'S0011'	3;...
        'S0012'	2;...
        'S0013'	2;...
        'S0014'	3;...
        'S0015'	3;...
        'S0016'	1;...
        'S0017'	2;...
        'S0018'	1;...
        'S0019'	1;...
        'S0020'	1;...
        'S0021'	3;...
        'S0022'	1;...
        'S0023'	3;...
        'S0024'	1;...
        'S0025'	2;...
        'S0026'	2;...
        'S0027'	1;...
        'S0028'	1;...
        'S0029'	2;...
        'S0030'	2;...
        'S0031'	2;...
        'S0032'	2;...
        'S0033'	2;...
        'S0034'	2;...
        'S0035'	2;...
        'S0036'	2;...
        'S0037'	1;...
        'S0038'	2;...
        'S0039'	2;...
        'S0040'	1;...
        'S0041'	1;...
        'S0042'	4;...
        'S0043'	4;...
        'S0044'	4;...
        'S0045'	4;...
        'S0046'	4;...
        'S0047'	4;...
        'S0048'	4;...
        'S0049'	4;...
        'S0050'	4;...
        'S0051'	4;...
        'S0052'	4;...
        'S0053'	2;...
        'S0054'	4;...
        'S0055'	4;...
        'S0056'	3;...
        'S0057'	3;...
        'S0058'	3;...
        'S0059'	3;...
        'S0060'	3 ...
        };
end
%% DECLARATION OF SCRIPT VARIABLES
workspace_file = fullfile(save_path,'dmri_comparison.mat');
if read_data == 1
    sbj = cell2mat(subject(:,2));

    protocol = {'ZOOMit_interp'; 'ZOOMit_notinterp' ;'RESOLVE'};
    name={'zoomit'; 'zoomit'; 'resolve'};

    cor_protocols = zeros(1,3);
    t_value = cor_protocols;
    p_value = cor_protocols;
    FA_gauss_gm = zeros(size(subject,1),3,size(protocol,1));
    FA_gauss_wm = zeros(size(subject,1),3,size(protocol,1));
    FAx_vec = 0:0.01:1;
    f1x_vec = 0:0.01:1;
    MDx_vec = 0:0.03:3;
    dx_vec = 0:0.04:4;
    FA_gm_mean = zeros(size(subject,1),size(protocol,1));
    FA_wm_mean = zeros(size(subject,1),size(protocol,1));
    FA_gm_var = zeros(size(subject,1),size(protocol,1));
    FA_wm_var = zeros(size(subject,1),size(protocol,1));
    FA_gm_median = zeros(size(subject,1),size(protocol,1));
    FA_wm_median = zeros(size(subject,1),size(protocol,1));
    FA_gm_mode = zeros(size(subject,1),size(protocol,1));
    FA_wm_mode = zeros(size(subject,1),size(protocol,1));
    FA_gm_skewness = zeros(size(subject,1),size(protocol,1));
    FA_wm_skewness = zeros(size(subject,1),size(protocol,1));
    FA_gm_kurtosis = zeros(size(subject,1),size(protocol,1));
    FA_wm_kurtosis = zeros(size(subject,1),size(protocol,1));
    FA_gm_kernel = zeros(size(FAx_vec,2),size(subject,1),size(protocol,1));
    FA_wm_kernel = zeros(size(FAx_vec,2),size(subject,1),size(protocol,1));
    FA_diff_wm_gm_median = zeros(size(subject,1),size(protocol,1));
    FA_diff_wm_gm_mean = zeros(size(subject,1),size(protocol,1));
    FA_diff_wm_gm_mode = zeros(size(subject,1),size(protocol,1));
    FA_gm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    FA_wm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    FA_gm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    FA_wm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    f1_gauss_gm = zeros(size(subject,1),3,size(protocol,1));
    f1_gauss_wm = zeros(size(subject,1),3,size(protocol,1));
    f1_gauss_gm_mean = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_mean = zeros(size(subject,1),size(protocol,1));
    f1_gauss_gm_var = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_var = zeros(size(subject,1),size(protocol,1));
    f1_gauss_gm_median = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_median = zeros(size(subject,1),size(protocol,1));
    f1_gauss_gm_mode = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_mode = zeros(size(subject,1),size(protocol,1));
    f1_gauss_gm_skewness = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_skewness = zeros(size(subject,1),size(protocol,1));
    f1_gauss_gm_kurtosis = zeros(size(subject,1),size(protocol,1));
    f1_gauss_wm_kurtosis = zeros(size(subject,1),size(protocol,1));
    f1_gm_kernel = zeros(size(f1x_vec,2),size(subject,1),size(protocol,1));
    f1_wm_kernel = zeros(size(f1x_vec,2),size(subject,1),size(protocol,1));
    f1_gauss_diff_wm_gm_median = zeros(size(subject,1),size(protocol,1));
    f1_gauss_diff_wm_gm_mean = zeros(size(subject,1),size(protocol,1));
    f1_gauss_diff_wm_gm_mode = zeros(size(subject,1),size(protocol,1));
    f1_gm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    f1_wm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    f1_gm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    f1_wm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    MD_gauss_gm = zeros(size(subject,1),3,size(protocol,1));
    MD_gauss_wm = zeros(size(subject,1),3,size(protocol,1));
    MD_gauss_gm_mean = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_mean = zeros(size(subject,1),size(protocol,1));
    MD_gauss_gm_var = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_var = zeros(size(subject,1),size(protocol,1));
    MD_gauss_gm_median = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_median = zeros(size(subject,1),size(protocol,1));
    MD_gauss_gm_mode = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_mode = zeros(size(subject,1),size(protocol,1));
    MD_gauss_gm_skewness = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_skewness = zeros(size(subject,1),size(protocol,1));
    MD_gauss_gm_kurtosis = zeros(size(subject,1),size(protocol,1));
    MD_gauss_wm_kurtosis = zeros(size(subject,1),size(protocol,1));
    MD_gm_kernel = zeros(size(MDx_vec,2),size(subject,1),size(protocol,1));
    MD_wm_kernel = zeros(size(MDx_vec,2),size(subject,1),size(protocol,1));
    MD_gauss_diff_wm_gm_median = zeros(size(subject,1),size(protocol,1));
    MD_gauss_diff_wm_gm_mean = zeros(size(subject,1),size(protocol,1));
    MD_gauss_diff_wm_gm_mode = zeros(size(subject,1),size(protocol,1));
    MD_gm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    MD_wm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    MD_gm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    MD_wm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    d_gauss_gm = zeros(size(subject,1),3,size(protocol,1));
    d_gauss_wm = zeros(size(subject,1),3,size(protocol,1));
    d_gauss_gm_mean = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_mean = zeros(size(subject,1),size(protocol,1));
    d_gauss_gm_var = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_var = zeros(size(subject,1),size(protocol,1));
    d_gauss_gm_median = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_median = zeros(size(subject,1),size(protocol,1));
    d_gauss_gm_mode = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_mode = zeros(size(subject,1),size(protocol,1));
    d_gauss_gm_skewness = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_skewness = zeros(size(subject,1),size(protocol,1));
    d_gauss_gm_kurtosis = zeros(size(subject,1),size(protocol,1));
    d_gauss_wm_kurtosis = zeros(size(subject,1),size(protocol,1));
    d_gm_kernel = zeros(size(dx_vec,2),size(subject,1),size(protocol,1));
    d_wm_kernel = zeros(size(dx_vec,2),size(subject,1),size(protocol,1));
    d_gauss_diff_wm_gm_median = zeros(size(subject,1),size(protocol,1));
    d_gauss_diff_wm_gm_mean = zeros(size(subject,1),size(protocol,1));
    d_gauss_diff_wm_gm_mode = zeros(size(subject,1),size(protocol,1));
    d_gm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    d_wm_pKolSmir = 3*ones(size(subject,1),size(protocol,1));
    d_gm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    d_wm_hKolSmir = 3*ones(size(subject,1),size(protocol,1));
    spl = 1;
    fieldcoef_mean = zeros(size(subject,1),size(protocol,1));
    fieldcoef_median = zeros(size(subject,1),size(protocol,1));
    fieldcoef_var = zeros(size(subject,1),size(protocol,1));
    fieldcoef3_mean = zeros(size(subject,1),size(protocol,1));
    fieldcoef3_median = zeros(size(subject,1),size(protocol,1));
    fieldcoef3_var = zeros(size(subject,1),size(protocol,1));
    fieldcoef56_mean = zeros(size(subject,1),size(protocol,1));
    fieldcoef56_median = zeros(size(subject,1),size(protocol,1));
    fieldcoef56_var = zeros(size(subject,1),size(protocol,1));
    
    if plot_results == 1
        h(1).fig = figure(1);
        set(h(1).fig, 'Position', [50, 50, 850, 1100]);
    end
%% DMRI DESCRIPTIVE STATISTICS METRIC EXTRACTION
    for sub = 1:size(subject,1)
        for ind = 1:size(protocol,1)
            if ind > 4
                dmri_preproc_folder = fullfile(data_folder,subject{sub,1},'Results',['Diff_preproc_' protocol{ind,1} '_moco']);
            else
                dmri_preproc_folder = fullfile(data_folder,subject{sub,1},'Results',['Diff_preproc_' protocol{ind,1}]);
            end

            anat_preproc_folder = fullfile(data_folder,subject{sub,1},'Results','Anat_Preproc');

            T2TRA_file=fullfile(anat_preproc_folder,'T2TRA_thr_bias_corr.nii');
            gunzip([T2TRA_file '.gz']);
            T2TRA_file_hdr = spm_vol(T2TRA_file);
            T2TRA = spm_read_vols(T2TRA_file_hdr);
            delete(T2TRA_file);

            label_file = fullfile(anat_preproc_folder,'T2TRA_thr_bias_corr_seg_labeled.nii');
            gunzip([label_file '.gz']);
            label_file_hdr = spm_vol(label_file);
            label = spm_read_vols(label_file_hdr);
            delete(label_file);

            gm_file = fullfile(anat_preproc_folder,'T2TRA_thr_bias_corr_gmseg.nii');
            gunzip([gm_file '.gz']);
            gm_file_hdr = spm_vol(gm_file);
            gm = spm_read_vols(gm_file_hdr);
            delete(gm_file);

            mask = zeros(size(T2TRA));
            mask(label >= 3 & label <= 6) = 1;

            mask_gm = zeros(size(T2TRA));
            mask_gm(label >= 3 & label <= 6 & gm == 1) = 1;

            mask_wm = mask - mask_gm;

            mask3 = zeros(size(T2TRA));
            mask3(label == 3) = 1;

            mask56 = zeros(size(T2TRA));
            mask56(label >= 5 & label <= 6) = 1;

            if isfolder(dmri_preproc_folder)                           
                FA_file=fullfile(dmri_preproc_folder,'diff_in_t2tra',['dti_' name{ind,1} '_FA.nii']);
                gunzip([FA_file '.gz']);
                FA_file_hdr = spm_vol(FA_file);
                FA = spm_read_vols(FA_file_hdr);
                delete(FA_file)

                MD_file=fullfile(dmri_preproc_folder,'diff_in_t2tra',['dti_' name{ind,1} '_MD.nii']);
                pom_file=[MD_file '.gz'];
                gunzip(pom_file);
                MD_file_hdr = spm_vol(MD_file);
                MD = spm_read_vols(MD_file_hdr);
                delete(MD_file)

                f1_file = fullfile(dmri_preproc_folder,'diff_in_t2tra','mean_f1samples.nii');
                gunzip([f1_file '.gz']);
                f1_file_hdr = spm_vol(f1_file);
                f1 = spm_read_vols(f1_file_hdr);
                delete(f1_file);

                f2_file = fullfile(dmri_preproc_folder,'diff_in_t2tra','mean_f2samples.nii');
                gunzip([f2_file '.gz']);
                f2_file_hdr = spm_vol(f2_file);
                f2 = spm_read_vols(f2_file_hdr);
                delete(f2_file);

                d_file = fullfile(dmri_preproc_folder,'diff_in_t2tra','mean_dsamples.nii');
                gunzip([d_file '.gz']);
                d_file_hdr = spm_vol(d_file);
                d = spm_read_vols(d_file_hdr);
                delete(d_file)

                fieldcoef_file = fullfile(dmri_preproc_folder,'diff_in_t2tra',['field_' name{ind,1} '_topup.nii']);
                gunzip([fieldcoef_file '.gz']);
                fieldcoef_file_hdr = spm_vol(fieldcoef_file);
                fieldcoef = spm_read_vols(fieldcoef_file_hdr);
                delete(fieldcoef_file);

                FA_vec = FA(mask==1);
                T2TRA_vec = T2TRA(mask==1);
                f1_vec = f1(mask==1);
                f2_vec = f2(mask==1);
                f2_vec(f2_vec<0.05) = [];
                MD_vec = MD(mask==1);
                MD_vec_norm = (MD_vec - min(MD_vec)) / (max(MD_vec)-min(MD_vec));
                d_vec = d(mask==1);
                d_vec_norm = (d_vec - min(d_vec)) / (max(d_vec)-min(d_vec));
                fieldcoef_vec = abs(fieldcoef(mask==1));
                fieldcoef_vec3 = abs(fieldcoef(mask3==1));
                fieldcoef_vec56 = abs(fieldcoef(mask56==1));

                cor = corrcoef(FA_vec,T2TRA_vec);
                cor_protocols(sub,ind) = cor(1,2);

                n=size(FA_vec,1);
                t_value(sub,ind) = cor_protocols(sub,ind) / sqrt( (1-cor_protocols(sub,ind)^2)/(n-2)  );
                p_value(sub,ind) = 2*(1-tcdf(abs(t_value(sub,ind)),n-1));

                pom = [FA_vec T2TRA_vec];
                pom = ( pom - min(pom(:)) ) / ( max(pom(:)) - min(pom(:)) );
                pom = round(2^16*pom);
                MI(sub,ind) = nmi(pom(:,1),pom(:,2));

                f2_prob(sub,ind) = ( size(f2_vec,1) / sum(mask(:)) ) * 100;

                pom = [f1_vec T2TRA_vec];
                pom = ( pom - min(pom(:)) ) / ( max(pom(:)) - min(pom(:)) );
                pom = round(2^16*pom);
                MIf1(sub,ind) = nmi(pom(:,1),pom(:,2));

                pom = [MD_vec_norm T2TRA_vec];
                pom = ( pom - min(pom(:)) ) / ( max(pom(:)) - min(pom(:)) );
                pom = round(2^16*pom);
                MIMD(sub,ind) = nmi(pom(:,1),pom(:,2));

                pom = [d_vec_norm T2TRA_vec];
                pom = ( pom - min(pom(:)) ) / ( max(pom(:)) - min(pom(:)) );
                pom = round(2^16*pom);
                MId(sub,ind) = nmi(pom(:,1),pom(:,2));

                FA_gm = FA(mask_gm==1);
                FA_wm = FA(mask_wm==1);
                FA_gm = FA_gm(FA_gm>0 & FA_gm<=1);
                FA_wm = FA_wm(FA_wm>0 & FA_wm<=1);

                f1_gm = f1(mask_gm==1);
                f1_wm = f1(mask_wm==1);
                f1_gm = f1_gm(f1_gm>0 & f1_gm<=1);
                f1_wm = f1_wm(f1_wm>0 & f1_wm<=1);

                MD_gm = MD(mask_gm==1)*10^3;
                MD_wm = MD(mask_wm==1)*10^3;
                MD_gm = MD_gm(MD_gm>0);
                MD_wm = MD_wm(MD_wm>0);

                d_gm = d(mask_gm==1)*10^3;
                d_wm = d(mask_wm==1)*10^3;
                d_gm = d_gm(d_gm>0);
                d_wm = d_wm(d_wm>0);

                FA_gm_mean(sub,ind) = mean(FA_gm);
                FA_wm_mean(sub,ind) = mean(FA_wm);
                FA_gm_var(sub,ind) = var(FA_gm);
                FA_wm_var(sub,ind) = var(FA_wm);
                FA_gm_median(sub,ind) = median(FA_gm);
                FA_wm_median(sub,ind) = median(FA_wm);
                FA_gm_skewness(sub,ind) = skewness(FA_gm);
                FA_wm_skewness(sub,ind) = skewness(FA_wm);
                FA_gm_kurtosis(sub,ind) = kurtosis(FA_gm);
                FA_wm_kurtosis(sub,ind) = kurtosis(FA_wm);
                FA_diff_wm_gm_median(sub,ind) = FA_wm_median(sub,ind) - FA_gm_median(sub,ind);
                FA_diff_wm_gm_mean(sub,ind) = FA_wm_mean(sub,ind) - FA_gm_mean(sub,ind);
                pd = fitdist(FA_gm,'Kernel');
                y_pdf = pdf(pd,FAx_vec);
                FA_gm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                FA_gm_mode(sub,ind) = FAx_vec(y_pdf==max(y_pdf));
                pd = fitdist(FA_wm,'Kernel');
                y_pdf = pdf(pd,FAx_vec);
                FA_wm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                FA_wm_mode(sub,ind) = FAx_vec(y_pdf==max(y_pdf));
                FA_diff_wm_gm_mode(sub,ind) = FA_wm_mode(sub,ind) - FA_gm_mode(sub,ind);
                [FA_gm_hKolSmir(sub,ind), FA_gm_pKolSmir(sub,ind)] = kstest(FA_gm);
                [FA_wm_hKolSmir(sub,ind), FA_wm_pKolSmir(sub,ind)] = kstest(FA_wm);

                f1_gauss_gm_mean(sub,ind) = mean(f1_gm);
                f1_gauss_wm_mean(sub,ind) = mean(f1_wm);
                f1_gauss_gm_var(sub,ind) = var(f1_gm);
                f1_gauss_wm_var(sub,ind) = var(f1_wm);
                f1_gauss_gm_median(sub,ind) = median(f1_gm);
                f1_gauss_wm_median(sub,ind) = median(f1_wm);
                f1_gauss_gm_skewness(sub,ind) = skewness(f1_gm);
                f1_gauss_wm_skewness(sub,ind) = skewness(f1_wm);
                f1_gauss_gm_kurtosis(sub,ind) = kurtosis(f1_gm);
                f1_gauss_wm_kurtosis(sub,ind) = kurtosis(f1_wm);
                f1_gauss_diff_wm_gm_median(sub,ind) = f1_gauss_wm_median(sub,ind) - f1_gauss_gm_median(sub,ind);
                f1_gauss_diff_wm_gm_mean(sub,ind) = f1_gauss_wm_mean(sub,ind) - f1_gauss_gm_mean(sub,ind);
                pd = fitdist(f1_gm,'Kernel');
                y_pdf = pdf(pd,f1x_vec);
                f1_gm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                f1_gauss_gm_mode(sub,ind) = f1x_vec(y_pdf==max(y_pdf));
                pd = fitdist(f1_wm,'Kernel');
                y_pdf = pdf(pd,f1x_vec);
                f1_wm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                f1_gauss_wm_mode(sub,ind) = f1x_vec(y_pdf==max(y_pdf));
                f1_gauss_diff_wm_gm_mode(sub,ind) = f1_gauss_wm_mode(sub,ind) - f1_gauss_gm_mode(sub,ind);
                [f1_gm_hKolSmir(sub,ind), f1_gm_pKolSmir(sub,ind)] = kstest(f1_gm);
                [f1_wm_hKolSmir(sub,ind), f1_wm_pKolSmir(sub,ind)] = kstest(f1_wm);

                MD_gauss_gm_mean(sub,ind) = mean(MD_gm);
                MD_gauss_wm_mean(sub,ind) = mean(MD_wm);
                MD_gauss_gm_var(sub,ind) = var(MD_gm);
                MD_gauss_wm_var(sub,ind) = var(MD_wm);
                MD_gauss_gm_median(sub,ind) = median(MD_gm);
                MD_gauss_wm_median(sub,ind) = median(MD_wm);
                MD_gauss_gm_skewness(sub,ind) = skewness(MD_gm);
                MD_gauss_wm_skewness(sub,ind) = skewness(MD_wm);
                MD_gauss_gm_kurtosis(sub,ind) = kurtosis(MD_gm);
                MD_gauss_wm_kurtosis(sub,ind) = kurtosis(MD_wm);
                MD_gauss_diff_wm_gm_median(sub,ind) = MD_gauss_wm_median(sub,ind) - MD_gauss_gm_median(sub,ind);
                MD_gauss_diff_wm_gm_mean(sub,ind) = MD_gauss_wm_mean(sub,ind) - MD_gauss_gm_mean(sub,ind);
                pd = fitdist(MD_gm,'Kernel');
                y_pdf = pdf(pd,MDx_vec);
                MD_gm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                MD_gauss_gm_mode(sub,ind) = MDx_vec(y_pdf==max(y_pdf));
                pd = fitdist(MD_wm,'Kernel');
                y_pdf = pdf(pd,MDx_vec);
                MD_wm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                MD_gauss_wm_mode(sub,ind) = MDx_vec(y_pdf==max(y_pdf));
                MD_gauss_diff_wm_gm_mode(sub,ind) = MD_gauss_wm_mode(sub,ind) - MD_gauss_gm_mode(sub,ind);
                [MD_gm_hKolSmir(sub,ind), MD_gm_pKolSmir(sub,ind)] = kstest(MD_gm);
                [MD_wm_hKolSmir(sub,ind), MD_wm_pKolSmir(sub,ind)] = kstest(MD_wm);

                d_gauss_gm_mean(sub,ind) = mean(d_gm);
                d_gauss_wm_mean(sub,ind) = mean(d_wm);
                d_gauss_gm_var(sub,ind) = var(d_gm);
                d_gauss_wm_var(sub,ind) = var(d_wm);
                d_gauss_gm_median(sub,ind) = median(d_gm);
                d_gauss_wm_median(sub,ind) = median(d_wm);
                d_gauss_gm_skewness(sub,ind) = skewness(d_gm);
                d_gauss_wm_skewness(sub,ind) = skewness(d_wm);
                d_gauss_gm_kurtosis(sub,ind) = kurtosis(d_gm);
                d_gauss_wm_kurtosis(sub,ind) = kurtosis(d_wm);
                d_gauss_diff_wm_gm_median(sub,ind) = d_gauss_wm_median(sub,ind) - d_gauss_gm_median(sub,ind);
                d_gauss_diff_wm_gm_mean(sub,ind) = d_gauss_wm_mean(sub,ind) - d_gauss_gm_mean(sub,ind);
                pd = fitdist(d_gm,'Kernel');
                y_pdf = pdf(pd,dx_vec);
                d_gm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                d_gauss_gm_mode(sub,ind) = dx_vec(y_pdf==max(y_pdf));
                pd = fitdist(d_wm,'Kernel');
                y_pdf = pdf(pd,dx_vec);
                d_wm_kernel(:,sub,ind) = y_pdf/sum(y_pdf);
                d_gauss_wm_mode(sub,ind) = dx_vec(y_pdf==max(y_pdf));
                d_gauss_diff_wm_gm_mode(sub,ind) = d_gauss_wm_mode(sub,ind) - d_gauss_gm_mode(sub,ind);
                [d_gm_hKolSmir(sub,ind), d_gm_pKolSmir(sub,ind)] = kstest(d_gm);
                [d_wm_hKolSmir(sub,ind), d_wm_pKolSmir(sub,ind)] = kstest(d_wm);

                fieldcoef_mean(sub,ind) = mean(fieldcoef_vec);
                fieldcoef_median(sub,ind) = median(fieldcoef_vec);
                if fieldcoef_median(sub,ind) == 0
                    fieldcoef_median(sub,ind) = 1e-9;
                end
                fieldcoef_var(sub,ind) = var(fieldcoef_vec);

                fieldcoef3_mean(sub,ind) = mean(fieldcoef_vec3);
                fieldcoef3_median(sub,ind) = median(fieldcoef_vec3);
                if fieldcoef3_median(sub,ind) == 0
                    fieldcoef3_median(sub,ind) = 1e-9;
                end
                fieldcoef3_var(sub,ind) = var(fieldcoef_vec3);

                fieldcoef56_mean(sub,ind) = mean(fieldcoef_vec56);
                fieldcoef56_median(sub,ind) = median(fieldcoef_vec56);
                if fieldcoef56_median(sub,ind) == 0
                    fieldcoef56_median(sub,ind) = 1e-9;
                end
                fieldcoef56_var(sub,ind) = var(fieldcoef_vec56);

                if plot_results == 1    
                    if sub<=8
                        subplot(8,3,spl)
                        H1 = histogram(FA_wm,35,'Normalization','probability');
                        hold on
                        H2 = histogram(FA_gm,35,'Normalization','probability');
                        FA_gm_pdf = FA_gm_kernel(:,sub,ind)';
                        FA_gm_pdf = (FA_gm_pdf/max(FA_gm_pdf)) * 0.965*max(H2.BinCounts/sum(H2.BinCounts));
                        FA_wm_pdf = FA_wm_kernel(:,sub,ind)';
                        FA_wm_pdf = (FA_wm_pdf/max(FA_wm_pdf)) * 0.965*max(H1.BinCounts/sum(H1.BinCounts));
                        stem(FA_wm_median(sub,ind),0.95*max(H1.BinCounts/sum(H1.BinCounts)),'b-.','LineWidth',8,'Marker','none')
                        stem(FA_gm_median(sub,ind),0.95*max(H2.BinCounts/sum(H2.BinCounts)),'r-.','LineWidth',8,'Marker','none')
                        H3 = plot(FAx_vec,FA_wm_pdf,'b-.','LineWidth',8);
                        H4 = plot(FAx_vec,FA_gm_pdf,'r-.','LineWidth',8);
                        hold off
                        xlim([0 1]);
                        ylim([0 1.05*max([FA_wm_pdf FA_gm_pdf])/0.965])
                        if ind == 1 && subject{sub,2}==1
                            ylabel('Control')
                        elseif ind == 1 && subject{sub,2}==2
                            ylabel('Low comp.')
                        elseif ind == 1 && subject{sub,2}==3
                            ylabel('Reproducibility')
                        elseif ind == 1 && subject{sub,2}==4
                            ylabel('Strong comp.')
                        end
                        if sub == 8 && ind == 1
                            xlabel('ZOOMit interp')
                        elseif sub == 8 && ind == 2
                            xlabel('ZOOMit notinterp')
                        elseif sub == 8 && ind == 3
                            xlabel('RESOLVE')
                        end
                        set(gca,'FontSize',12);
                        if sub == 8 && ind == 3
                            legend([H1 H2 H3 H4],{'WM','GM','WM smooth-dist.','GM smooth-dist.'});
                            set(h(1).fig.Children(1,1),'units','pixel')
                            set(h(1).fig.Children(1,1),'FontSize',12)
                            set(h(1).fig.Children(1,1),'position',[720 1000 110 70])

                            set(gcf,'NextPlot','add');
                            axes;
                            nadpis = title({'Estimated FA Normal distributions with medians from C3-C6 WM and GM';'x-axis - FA values; y-axis - probability'},'FontSize',14);
                            set(gca,'Visible','off');
                            set(nadpis,'Visible','on');

                            pause(1)
                            print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_WM_GM_distributions']), '-dpng', '-r300')
                        end
                        pause(0.1)

                        spl = spl + 1;
                    end
                end
            end
            disp(['Subject ' num2str(sub) ' session ' num2str(ind) ' done.'])
        end
    end
%% IINSERTING RESULTS INTO DMRI PROTOCOL-SPECIFIC TABLES
    table_ZOOMit_Int(:,1) = FA_wm_mean(:,1);
    table_ZOOMit_Int(:,2) = FA_wm_median(:,1);
    table_ZOOMit_Int(:,3) = sqrt(FA_wm_var(:,1));
    table_ZOOMit_Int(:,4) = FA_wm_mode(:,1);
    table_ZOOMit_Int(:,5) = FA_gm_mean(:,1);
    table_ZOOMit_Int(:,6) = FA_gm_median(:,1);
    table_ZOOMit_Int(:,7) = sqrt(FA_gm_var(:,1));
    table_ZOOMit_Int(:,8) = FA_gm_mode(:,1);
    table_ZOOMit_Int(:,9) = FA_diff_wm_gm_mean(:,1);
    table_ZOOMit_Int(:,10) = FA_diff_wm_gm_median(:,1);
    table_ZOOMit_Int(:,11) = MD_gauss_wm_mean(:,1);
    table_ZOOMit_Int(:,12) = MD_gauss_wm_median(:,1);
    table_ZOOMit_Int(:,13) = sqrt(MD_gauss_wm_var(:,1));
    table_ZOOMit_Int(:,14) = MD_gauss_wm_mode(:,1);
    table_ZOOMit_Int(:,15) = MD_gauss_gm_mean(:,1);
    table_ZOOMit_Int(:,16) = MD_gauss_gm_median(:,1);
    table_ZOOMit_Int(:,17) = sqrt(MD_gauss_gm_var(:,1));
    table_ZOOMit_Int(:,18) = MD_gauss_gm_mode(:,1);
    table_ZOOMit_Int(:,19) = MD_gauss_diff_wm_gm_mean(:,1);
    table_ZOOMit_Int(:,20) = MD_gauss_diff_wm_gm_median(:,1);
    table_ZOOMit_Int(:,21) = f1_gauss_wm_mean(:,1);
    table_ZOOMit_Int(:,22) = f1_gauss_wm_median(:,1);
    table_ZOOMit_Int(:,23) = sqrt(f1_gauss_wm_var(:,1));
    table_ZOOMit_Int(:,24) = f1_gauss_wm_mode(:,1);
    table_ZOOMit_Int(:,25) = f1_gauss_gm_mean(:,1);
    table_ZOOMit_Int(:,26) = f1_gauss_gm_median(:,1);
    table_ZOOMit_Int(:,27) = sqrt(f1_gauss_gm_var(:,1));
    table_ZOOMit_Int(:,28) = f1_gauss_gm_mode(:,1);
    table_ZOOMit_Int(:,29) = f1_gauss_diff_wm_gm_mean(:,1);
    table_ZOOMit_Int(:,30) = f1_gauss_diff_wm_gm_median(:,1);
    table_ZOOMit_Int(:,31) = d_gauss_wm_mean(:,1);
    table_ZOOMit_Int(:,32) = d_gauss_wm_median(:,1);
    table_ZOOMit_Int(:,33) = sqrt(d_gauss_wm_var(:,1));
    table_ZOOMit_Int(:,34) = d_gauss_wm_mode(:,1);
    table_ZOOMit_Int(:,35) = d_gauss_gm_mean(:,1);
    table_ZOOMit_Int(:,36) = d_gauss_gm_median(:,1);
    table_ZOOMit_Int(:,37) = sqrt(d_gauss_gm_var(:,1));
    table_ZOOMit_Int(:,38) = d_gauss_gm_mode(:,1);
    table_ZOOMit_Int(:,39) = d_gauss_diff_wm_gm_mean(:,1);
    table_ZOOMit_Int(:,40) = d_gauss_diff_wm_gm_median(:,1);
    table_ZOOMit_Int(:,41) = FA_wm_skewness(:,1);
    table_ZOOMit_Int(:,42) = MD_gauss_wm_skewness(:,1);
    table_ZOOMit_Int(:,43) = f1_gauss_wm_skewness(:,1);
    table_ZOOMit_Int(:,44) = d_gauss_wm_skewness(:,1);
    table_ZOOMit_Int(:,45) = FA_gm_skewness(:,1);
    table_ZOOMit_Int(:,46) = MD_gauss_gm_skewness(:,1);
    table_ZOOMit_Int(:,47) = f1_gauss_gm_skewness(:,1);
    table_ZOOMit_Int(:,48) = d_gauss_gm_skewness(:,1);
    table_ZOOMit_Int(:,49) = FA_wm_kurtosis(:,1);
    table_ZOOMit_Int(:,50) = MD_gauss_wm_kurtosis(:,1);
    table_ZOOMit_Int(:,51) = f1_gauss_wm_kurtosis(:,1);
    table_ZOOMit_Int(:,52) = d_gauss_wm_kurtosis(:,1);
    table_ZOOMit_Int(:,53) = FA_gm_kurtosis(:,1);
    table_ZOOMit_Int(:,54) = MD_gauss_gm_kurtosis(:,1);
    table_ZOOMit_Int(:,55) = f1_gauss_gm_kurtosis(:,1);
    table_ZOOMit_Int(:,56) = d_gauss_gm_kurtosis(:,1);
    table_ZOOMit_Int(:,57) = fieldcoef_mean(:,1);
    table_ZOOMit_Int(:,58) = fieldcoef3_mean(:,1);
    table_ZOOMit_Int(:,59) = fieldcoef56_mean(:,1);
    table_ZOOMit_Int(:,60) = fieldcoef_median(:,1);
    table_ZOOMit_Int(:,61) = fieldcoef3_median(:,1);
    table_ZOOMit_Int(:,62) = fieldcoef56_median(:,1);
    table_ZOOMit_Int(:,63) = sqrt(fieldcoef_var(:,1));
    table_ZOOMit_Int(:,64) = sqrt(fieldcoef3_var(:,1));
    table_ZOOMit_Int(:,65) = sqrt(fieldcoef56_var(:,1));
    table_ZOOMit_Int(:,66) = FA_diff_wm_gm_mode(:,1);
    table_ZOOMit_Int(:,67) = f1_gauss_diff_wm_gm_mode(:,1);
    table_ZOOMit_Int(:,68) = MD_gauss_diff_wm_gm_mode(:,1);
    table_ZOOMit_Int(:,69) = d_gauss_diff_wm_gm_mode(:,1);

    table_ZOOMit_NotInt(:,1) = FA_wm_mean(:,2);
    table_ZOOMit_NotInt(:,2) = FA_wm_median(:,2);
    table_ZOOMit_NotInt(:,3) = sqrt(FA_wm_var(:,2));
    table_ZOOMit_NotInt(:,4) = FA_wm_mode(:,2);
    table_ZOOMit_NotInt(:,5) = FA_gm_mean(:,2);
    table_ZOOMit_NotInt(:,6) = FA_gm_median(:,2);
    table_ZOOMit_NotInt(:,7) = sqrt(FA_gm_var(:,2));
    table_ZOOMit_NotInt(:,8) = FA_gm_mode(:,2);
    table_ZOOMit_NotInt(:,9) = FA_diff_wm_gm_mean(:,2);
    table_ZOOMit_NotInt(:,10) = FA_diff_wm_gm_median(:,2);
    table_ZOOMit_NotInt(:,11) = MD_gauss_wm_mean(:,2);
    table_ZOOMit_NotInt(:,12) = MD_gauss_wm_median(:,2);
    table_ZOOMit_NotInt(:,13) = sqrt(MD_gauss_wm_var(:,2));
    table_ZOOMit_NotInt(:,14) = MD_gauss_wm_mode(:,2);
    table_ZOOMit_NotInt(:,15) = MD_gauss_gm_mean(:,2);
    table_ZOOMit_NotInt(:,16) = MD_gauss_gm_median(:,2);
    table_ZOOMit_NotInt(:,17) = sqrt(MD_gauss_gm_var(:,2));
    table_ZOOMit_NotInt(:,18) = MD_gauss_gm_mode(:,2);
    table_ZOOMit_NotInt(:,19) = MD_gauss_diff_wm_gm_mean(:,2);
    table_ZOOMit_NotInt(:,20) = MD_gauss_diff_wm_gm_median(:,2);
    table_ZOOMit_NotInt(:,21) = f1_gauss_wm_mean(:,2);
    table_ZOOMit_NotInt(:,22) = f1_gauss_wm_median(:,2);
    table_ZOOMit_NotInt(:,23) = sqrt(f1_gauss_wm_var(:,2));
    table_ZOOMit_NotInt(:,24) = f1_gauss_wm_mode(:,2);
    table_ZOOMit_NotInt(:,25) = f1_gauss_gm_mean(:,2);
    table_ZOOMit_NotInt(:,26) = f1_gauss_gm_median(:,2);
    table_ZOOMit_NotInt(:,27) = sqrt(f1_gauss_gm_var(:,2));
    table_ZOOMit_NotInt(:,28) = f1_gauss_gm_mode(:,2);
    table_ZOOMit_NotInt(:,29) = f1_gauss_diff_wm_gm_mean(:,2);
    table_ZOOMit_NotInt(:,30) = f1_gauss_diff_wm_gm_median(:,2);
    table_ZOOMit_NotInt(:,31) = d_gauss_wm_mean(:,2);
    table_ZOOMit_NotInt(:,32) = d_gauss_wm_median(:,2);
    table_ZOOMit_NotInt(:,33) = sqrt(d_gauss_wm_var(:,2));
    table_ZOOMit_NotInt(:,34) = d_gauss_wm_mode(:,2);
    table_ZOOMit_NotInt(:,35) = d_gauss_gm_mean(:,2);
    table_ZOOMit_NotInt(:,36) = d_gauss_gm_median(:,2);
    table_ZOOMit_NotInt(:,37) = sqrt(d_gauss_gm_var(:,2));
    table_ZOOMit_NotInt(:,38) = d_gauss_gm_mode(:,2);
    table_ZOOMit_NotInt(:,39) = d_gauss_diff_wm_gm_mean(:,2);
    table_ZOOMit_NotInt(:,40) = d_gauss_diff_wm_gm_median(:,2);
    table_ZOOMit_NotInt(:,41) = FA_wm_skewness(:,2);
    table_ZOOMit_NotInt(:,42) = MD_gauss_wm_skewness(:,2);
    table_ZOOMit_NotInt(:,43) = f1_gauss_wm_skewness(:,2);
    table_ZOOMit_NotInt(:,44) = d_gauss_wm_skewness(:,2);
    table_ZOOMit_NotInt(:,45) = FA_gm_skewness(:,2);
    table_ZOOMit_NotInt(:,46) = MD_gauss_gm_skewness(:,2);
    table_ZOOMit_NotInt(:,47) = f1_gauss_gm_skewness(:,2);
    table_ZOOMit_NotInt(:,48) = d_gauss_gm_skewness(:,2);
    table_ZOOMit_NotInt(:,49) = FA_wm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,50) = MD_gauss_wm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,51) = f1_gauss_wm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,52) = d_gauss_wm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,53) = FA_gm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,54) = MD_gauss_gm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,55) = f1_gauss_gm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,56) = d_gauss_gm_kurtosis(:,2);
    table_ZOOMit_NotInt(:,57) = fieldcoef_mean(:,2);
    table_ZOOMit_NotInt(:,58) = fieldcoef3_mean(:,2);
    table_ZOOMit_NotInt(:,59) = fieldcoef56_mean(:,2);
    table_ZOOMit_NotInt(:,60) = fieldcoef_median(:,2);
    table_ZOOMit_NotInt(:,61) = fieldcoef3_median(:,2);
    table_ZOOMit_NotInt(:,62) = fieldcoef56_median(:,2);
    table_ZOOMit_NotInt(:,63) = sqrt(fieldcoef_var(:,2));
    table_ZOOMit_NotInt(:,64) = sqrt(fieldcoef3_var(:,2));
    table_ZOOMit_NotInt(:,65) = sqrt(fieldcoef56_var(:,2));
    table_ZOOMit_NotInt(:,66) = FA_diff_wm_gm_mode(:,2);
    table_ZOOMit_NotInt(:,67) = f1_gauss_diff_wm_gm_mode(:,2);
    table_ZOOMit_NotInt(:,68) = MD_gauss_diff_wm_gm_mode(:,2);
    table_ZOOMit_NotInt(:,69) = d_gauss_diff_wm_gm_mode(:,2);

    table_RESOLVE(:,1) = FA_wm_mean(:,3);
    table_RESOLVE(:,2) = FA_wm_median(:,3);
    table_RESOLVE(:,3) = sqrt(FA_wm_var(:,3));
    table_RESOLVE(:,4) = FA_wm_mode(:,3);
    table_RESOLVE(:,5) = FA_gm_mean(:,3);
    table_RESOLVE(:,6) = FA_gm_median(:,3);
    table_RESOLVE(:,7) = sqrt(FA_gm_var(:,3));
    table_RESOLVE(:,8) = FA_gm_mode(:,3);
    table_RESOLVE(:,9) = FA_diff_wm_gm_mean(:,3);
    table_RESOLVE(:,10) = FA_diff_wm_gm_median(:,3);
    table_RESOLVE(:,11) = MD_gauss_wm_mean(:,3);
    table_RESOLVE(:,12) = MD_gauss_wm_median(:,3);
    table_RESOLVE(:,13) = sqrt(MD_gauss_wm_var(:,3));
    table_RESOLVE(:,14) = MD_gauss_wm_mode(:,3);
    table_RESOLVE(:,15) = MD_gauss_gm_mean(:,3);
    table_RESOLVE(:,16) = MD_gauss_gm_median(:,3);
    table_RESOLVE(:,17) = sqrt(MD_gauss_gm_var(:,3));
    table_RESOLVE(:,18) = MD_gauss_gm_mode(:,3);
    table_RESOLVE(:,19) = MD_gauss_diff_wm_gm_mean(:,3);
    table_RESOLVE(:,20) = MD_gauss_diff_wm_gm_median(:,3);
    table_RESOLVE(:,21) = f1_gauss_wm_mean(:,3);
    table_RESOLVE(:,22) = f1_gauss_wm_median(:,3);
    table_RESOLVE(:,23) = sqrt(f1_gauss_wm_var(:,3));
    table_RESOLVE(:,24) = f1_gauss_wm_mode(:,3);
    table_RESOLVE(:,25) = f1_gauss_gm_mean(:,3);
    table_RESOLVE(:,26) = f1_gauss_gm_median(:,3);
    table_RESOLVE(:,27) = sqrt(f1_gauss_gm_var(:,3));
    table_RESOLVE(:,28) = f1_gauss_gm_mode(:,3);
    table_RESOLVE(:,29) = f1_gauss_diff_wm_gm_mean(:,3);
    table_RESOLVE(:,30) = f1_gauss_diff_wm_gm_median(:,3);
    table_RESOLVE(:,31) = d_gauss_wm_mean(:,3);
    table_RESOLVE(:,32) = d_gauss_wm_median(:,3);
    table_RESOLVE(:,33) = sqrt(d_gauss_wm_var(:,3));
    table_RESOLVE(:,34) = d_gauss_wm_mode(:,3);
    table_RESOLVE(:,35) = d_gauss_gm_mean(:,3);
    table_RESOLVE(:,36) = d_gauss_gm_median(:,3);
    table_RESOLVE(:,37) = sqrt(d_gauss_gm_var(:,3));
    table_RESOLVE(:,38) = d_gauss_gm_mode(:,3);
    table_RESOLVE(:,39) = d_gauss_diff_wm_gm_mean(:,3);
    table_RESOLVE(:,40) = d_gauss_diff_wm_gm_median(:,3);
    table_RESOLVE(:,41) = FA_wm_skewness(:,3);
    table_RESOLVE(:,42) = MD_gauss_wm_skewness(:,3);
    table_RESOLVE(:,43) = f1_gauss_wm_skewness(:,3);
    table_RESOLVE(:,44) = d_gauss_wm_skewness(:,3);
    table_RESOLVE(:,45) = FA_gm_skewness(:,3);
    table_RESOLVE(:,46) = MD_gauss_gm_skewness(:,3);
    table_RESOLVE(:,47) = f1_gauss_gm_skewness(:,3);
    table_RESOLVE(:,48) = d_gauss_gm_skewness(:,3);
    table_RESOLVE(:,49) = FA_wm_kurtosis(:,3);
    table_RESOLVE(:,50) = MD_gauss_wm_kurtosis(:,3);
    table_RESOLVE(:,51) = f1_gauss_wm_kurtosis(:,3);
    table_RESOLVE(:,52) = d_gauss_wm_kurtosis(:,3);
    table_RESOLVE(:,53) = FA_gm_kurtosis(:,3);
    table_RESOLVE(:,54) = MD_gauss_gm_kurtosis(:,3);
    table_RESOLVE(:,55) = f1_gauss_gm_kurtosis(:,3);
    table_RESOLVE(:,56) = d_gauss_gm_kurtosis(:,3);
    table_RESOLVE(:,57) = fieldcoef_mean(:,3);
    table_RESOLVE(:,58) = fieldcoef3_mean(:,3);
    table_RESOLVE(:,59) = fieldcoef56_mean(:,3);
    table_RESOLVE(:,60) = fieldcoef_median(:,3);
    table_RESOLVE(:,61) = fieldcoef3_median(:,3);
    table_RESOLVE(:,62) = fieldcoef56_median(:,3);
    table_RESOLVE(:,63) = sqrt(fieldcoef_var(:,3));
    table_RESOLVE(:,64) = sqrt(fieldcoef3_var(:,3));
    table_RESOLVE(:,65) = sqrt(fieldcoef56_var(:,3));
    table_RESOLVE(:,66) = FA_diff_wm_gm_mode(:,3);
    table_RESOLVE(:,67) = f1_gauss_diff_wm_gm_mode(:,3);
    table_RESOLVE(:,68) = MD_gauss_diff_wm_gm_mode(:,3);
    table_RESOLVE(:,69) = d_gauss_diff_wm_gm_mode(:,3);
%% MEAN, MEDIAND AND STD METRIC GROUP STATISTICS
    for category = 1:3
        par_mat = table_ZOOMit_Int(sbj==category,:);
        par_mat(sum(par_mat,2)==0,:) = [];
        table_ZOOMit_Int_mean(category,:) = mean(par_mat);
        table_ZOOMit_Int_median(category,:) = median(par_mat);
        table_ZOOMit_Int_std(category,:) = std(par_mat);

        par_mat = table_ZOOMit_NotInt(sbj==category,:);
        par_mat(sum(par_mat,2)==0,:) = [];
        table_ZOOMit_NotInt_mean(category,:) = mean(par_mat);
        table_ZOOMit_NotInt_median(category,:) = median(par_mat);
        table_ZOOMit_NotInt_std(category,:) = std(par_mat);

        par_mat = table_RESOLVE(sbj==category,:);
        par_mat(sum(par_mat,2)==0,:) = [];
        table_RESOLVE_mean(category,:) = mean(par_mat);
        table_RESOLVE_median(category,:) = median(par_mat);
        table_RESOLVE_std(category,:) = std(par_mat);
    end
%% DECLARING OF EMPTY VARIABLES USED FOR SCATTERPLOT VISUALIZATIONS
    cor_boxplot = [];
    cor_grp = [];
    cor_mean = zeros(1,size(cor_protocols,2));
    cor_median = zeros(1,size(cor_protocols,2));
    t_value_boxplot = [];
    t_value_grp = [];
    t_value_mean = zeros(1,size(cor_protocols,2));
    t_value_median = zeros(1,size(cor_protocols,2));
    MI_boxplot = [];
    MI_grp = [];
    MI_mean = zeros(1,size(cor_protocols,2));
    MI_median = zeros(1,size(cor_protocols,2));
    subject_grp = [];
    f2_prob_boxplot = [];
    f2_prob_grp = [];
    f2_prob_median = zeros(1,size(cor_protocols,2));
    f2_prob_mean = zeros(1,size(cor_protocols,2));
    MIf1_boxplot = [];
    MIf1_grp = [];
    MIf1_mean = zeros(1,size(cor_protocols,2));
    MIf1_median = zeros(1,size(cor_protocols,2));
    MIMD_boxplot = [];
    MIMD_grp = [];
    MIMD_mean = zeros(1,size(cor_protocols,2));
    MIMD_median = zeros(1,size(cor_protocols,2));
    MId_boxplot = [];
    MId_grp = [];
    MId_mean = zeros(1,size(cor_protocols,2));
    MId_median = zeros(1,size(cor_protocols,2));
    FAWM_boxplot = [];
    FAWM_grp = [];
    FAWM_mean = zeros(1,size(cor_protocols,2));
    FAWM_median = zeros(1,size(cor_protocols,2));
    FAGM_boxplot = [];
    FAGM_grp = [];
    FAGM_mean = zeros(1,size(cor_protocols,2));
    FAGM_median = zeros(1,size(cor_protocols,2));
    FAWMGMdiff_boxplot = [];
    FAWMGMdiff_grp = [];
    FAWMGMdiff_mean = zeros(1,size(cor_protocols,2));
    FAWMGMdiff_median = zeros(1,size(cor_protocols,2));
    FAWMGMdiffmean_boxplot = [];
    FAWMGMdiffmean_grp = [];
    FAWMGMdiffmean_mean = zeros(1,size(cor_protocols,2));
    FAWMGMdiffmean_median = zeros(1,size(cor_protocols,2));
    FAWMGMdiffmode_boxplot = [];
    FAWMGMdiffmode_grp = [];
    FAWMGMdiffmode_mean = zeros(1,size(cor_protocols,2));
    FAWMGMdiffmode_median = zeros(1,size(cor_protocols,2));
    FAwmstd_boxplot = [];
    FAwmstd_grp = [];
    FAwmstd_mean = zeros(1,size(cor_protocols,2));
    FAwmstd_median = zeros(1,size(cor_protocols,2));
    FAgmstd_boxplot = [];
    FAgmstd_grp = [];
    FAgmstd_mean = zeros(1,size(cor_protocols,2));
    FAgmstd_median = zeros(1,size(cor_protocols,2));
    FAWMmean_boxplot = [];
    FAWMmean_grp = [];
    FAWMmean_mean = zeros(1,size(cor_protocols,2));
    FAWMmean_median = zeros(1,size(cor_protocols,2));
    FAGMmean_boxplot = [];
    FAGMmean_grp = [];
    FAGMmean_mean = zeros(1,size(cor_protocols,2));
    FAGMmean_median = zeros(1,size(cor_protocols,2));
    FAWMmode_boxplot = [];
    FAWMmode_grp = [];
    FAWMmode_mean = zeros(1,size(cor_protocols,2));
    FAWMmode_median = zeros(1,size(cor_protocols,2));
    FAGMmode_boxplot = [];
    FAGMmode_grp = [];
    FAGMmode_mean = zeros(1,size(cor_protocols,2));
    FAGMmode_median = zeros(1,size(cor_protocols,2));
    FAWMskewness_boxplot = [];
    FAWMskewness_grp = [];
    FAWMskewness_mean = zeros(1,size(cor_protocols,2));
    FAWMskewness_median = zeros(1,size(cor_protocols,2));
    FAGMskewness_boxplot = [];
    FAGMskewness_grp = [];
    FAGMskewness_mean = zeros(1,size(cor_protocols,2));
    FAGMskewness_median = zeros(1,size(cor_protocols,2));
    FAWMkurtosis_boxplot = [];
    FAWMkurtosis_grp = [];
    FAWMkurtosis_mean = zeros(1,size(cor_protocols,2));
    FAWMkurtosis_median = zeros(1,size(cor_protocols,2));
    FAGMkurtosis_boxplot = [];
    FAGMkurtosis_grp = [];
    FAGMkurtosis_mean = zeros(1,size(cor_protocols,2));
    FAGMkurtosis_median = zeros(1,size(cor_protocols,2));
    FAgmpKolSmir_boxplot = [];
    FAgmpKolSmir_grp = [];
    FAgmpKolSmir_mean = zeros(1,size(cor_protocols,2));
    FAgmpKolSmir_median = zeros(1,size(cor_protocols,2));
    FAwmpKolSmir_boxplot = [];
    FAwmpKolSmir_grp = [];
    FAwmpKolSmir_mean = zeros(1,size(cor_protocols,2));
    FAwmpKolSmir_median = zeros(1,size(cor_protocols,2));
    MDWM_boxplot = [];
    MDWM_grp = [];
    MDWM_mean = zeros(1,size(cor_protocols,2));
    MDWM_median = zeros(1,size(cor_protocols,2));
    MDGM_boxplot = [];
    MDGM_grp = [];
    MDGM_mean = zeros(1,size(cor_protocols,2));
    MDGM_median = zeros(1,size(cor_protocols,2));
    MDWMmean_boxplot = [];
    MDWMmean_grp = [];
    MDWMmean_mean = zeros(1,size(cor_protocols,2));
    MDWMmean_median = zeros(1,size(cor_protocols,2));
    MDGMmean_boxplot = [];
    MDGMmean_grp = [];
    MDGMmean_mean = zeros(1,size(cor_protocols,2));
    MDGMmean_median = zeros(1,size(cor_protocols,2));
    MDWMmode_boxplot = [];
    MDWMmode_grp = [];
    MDWMmode_mean = zeros(1,size(cor_protocols,2));
    MDWMmode_median = zeros(1,size(cor_protocols,2));
    MDGMmode_boxplot = [];
    MDGMmode_grp = [];
    MDGMmode_mean = zeros(1,size(cor_protocols,2));
    MDGMmode_median = zeros(1,size(cor_protocols,2));
    MDWMskewness_boxplot = [];
    MDWMskewness_grp = [];
    MDWMskewness_mean = zeros(1,size(cor_protocols,2));
    MDWMskewness_median = zeros(1,size(cor_protocols,2));
    MDGMskewness_boxplot = [];
    MDGMskewness_grp = [];
    MDGMskewness_mean = zeros(1,size(cor_protocols,2));
    MDGMskewness_median = zeros(1,size(cor_protocols,2));
    MDWMkurtosis_boxplot = [];
    MDWMkurtosis_grp = [];
    MDWMkurtosis_mean = zeros(1,size(cor_protocols,2));
    MDWMkurtosis_median = zeros(1,size(cor_protocols,2));
    MDGMkurtosis_boxplot = [];
    MDGMkurtosis_grp = [];
    MDGMkurtosis_mean = zeros(1,size(cor_protocols,2));
    MDGMkurtosis_median = zeros(1,size(cor_protocols,2));
    MDWMGMdiffmean_boxplot = [];
    MDWMGMdiffmean_grp = [];
    MDWMGMdiffmean_mean = zeros(1,size(cor_protocols,2));
    MDWMGMdiffmean_median = zeros(1,size(cor_protocols,2));
    MDWMGMdiffmode_boxplot = [];
    MDWMGMdiffmode_grp = [];
    MDWMGMdiffmode_mean = zeros(1,size(cor_protocols,2));
    MDWMGMdiffmode_median = zeros(1,size(cor_protocols,2));
    MDWMGMdiff_boxplot = [];
    MDWMGMdiff_grp = [];
    MDWMGMdiff_mean = zeros(1,size(cor_protocols,2));
    MDWMGMdiff_median = zeros(1,size(cor_protocols,2));
    MDwmstd_boxplot = [];
    MDwmstd_grp = [];
    MDwmstd_mean = zeros(1,size(cor_protocols,2));
    MDwmstd_median = zeros(1,size(cor_protocols,2));
    MDgmstd_boxplot = [];
    MDgmstd_grp = [];
    MDgmstd_mean = zeros(1,size(cor_protocols,2));
    MDgmstd_median = zeros(1,size(cor_protocols,2));
    MDgmpKolSmir_boxplot = [];
    MDgmpKolSmir_grp = [];
    MDgmpKolSmir_mean = zeros(1,size(cor_protocols,2));
    MDgmpKolSmir_median = zeros(1,size(cor_protocols,2));
    MDwmpKolSmir_boxplot = [];
    MDwmpKolSmir_grp = [];
    MDwmpKolSmir_mean = zeros(1,size(cor_protocols,2));
    MDwmpKolSmir_median = zeros(1,size(cor_protocols,2));
    f1WM_boxplot = [];
    f1WM_grp = [];
    f1WM_mean = zeros(1,size(cor_protocols,2));
    f1WM_median = zeros(1,size(cor_protocols,2));
    f1GM_boxplot = [];
    f1GM_grp = [];
    f1GM_mean = zeros(1,size(cor_protocols,2));
    f1GM_median = zeros(1,size(cor_protocols,2));
    f1WMGMdiff_boxplot = [];
    f1WMGMdiff_grp = [];
    f1WMGMdiff_mean = zeros(1,size(cor_protocols,2));
    f1WMGMdiff_median = zeros(1,size(cor_protocols,2));
    f1WMGMdiffmean_boxplot = [];
    f1WMGMdiffmean_grp = [];
    f1WMGMdiffmean_mean = zeros(1,size(cor_protocols,2));
    f1WMGMdiffmean_median = zeros(1,size(cor_protocols,2));
    f1WMGMdiffmode_boxplot = [];
    f1WMGMdiffmode_grp = [];
    f1WMGMdiffmode_mean = zeros(1,size(cor_protocols,2));
    f1WMGMdiffmode_median = zeros(1,size(cor_protocols,2));
    f1wmstd_boxplot = [];
    f1wmstd_grp = [];
    f1wmstd_mean = zeros(1,size(cor_protocols,2));
    f1wmstd_median = zeros(1,size(cor_protocols,2));
    f1gmstd_boxplot = [];
    f1gmstd_grp = [];
    f1gmstd_mean = zeros(1,size(cor_protocols,2));
    f1gmstd_median = zeros(1,size(cor_protocols,2));
    f1WMmean_boxplot = [];
    f1WMmean_grp = [];
    f1WMmean_mean = zeros(1,size(cor_protocols,2));
    f1WMmean_median = zeros(1,size(cor_protocols,2));
    f1GMmean_boxplot = [];
    f1GMmean_grp = [];
    f1GMmean_mean = zeros(1,size(cor_protocols,2));
    f1GMmean_median = zeros(1,size(cor_protocols,2));
    f1WMmode_boxplot = [];
    f1WMmode_grp = [];
    f1WMmode_mean = zeros(1,size(cor_protocols,2));
    f1WMmode_median = zeros(1,size(cor_protocols,2));
    f1GMmode_boxplot = [];
    f1GMmode_grp = [];
    f1GMmode_mean = zeros(1,size(cor_protocols,2));
    f1GMmode_median = zeros(1,size(cor_protocols,2));
    f1WMskewness_boxplot = [];
    f1WMskewness_grp = [];
    f1WMskewness_mean = zeros(1,size(cor_protocols,2));
    f1WMskewness_median = zeros(1,size(cor_protocols,2));
    f1GMskewness_boxplot = [];
    f1GMskewness_grp = [];
    f1GMskewness_mean = zeros(1,size(cor_protocols,2));
    f1GMskewness_median = zeros(1,size(cor_protocols,2));
    f1WMkurtosis_boxplot = [];
    f1WMkurtosis_grp = [];
    f1WMkurtosis_mean = zeros(1,size(cor_protocols,2));
    f1WMkurtosis_median = zeros(1,size(cor_protocols,2));
    f1GMkurtosis_boxplot = [];
    f1GMkurtosis_grp = [];
    f1GMkurtosis_mean = zeros(1,size(cor_protocols,2));
    f1GMkurtosis_median = zeros(1,size(cor_protocols,2));
    dWM_boxplot = [];
    dWM_grp = [];
    dWM_mean = zeros(1,size(cor_protocols,2));
    dWM_median = zeros(1,size(cor_protocols,2));
    dGM_boxplot = [];
    dGM_grp = [];
    dGM_mean = zeros(1,size(cor_protocols,2));
    dGM_median = zeros(1,size(cor_protocols,2));
    dWMGMdiff_boxplot = [];
    dWMGMdiff_grp = [];
    dWMGMdiff_mean = zeros(1,size(cor_protocols,2));
    dWMGMdiff_median = zeros(1,size(cor_protocols,2));
    dWMmean_boxplot = [];
    dWMmean_grp = [];
    dWMmean_mean = zeros(1,size(cor_protocols,2));
    dWMmean_median = zeros(1,size(cor_protocols,2));
    dGMmean_boxplot = [];
    dGMmean_grp = [];
    dGMmean_mean = zeros(1,size(cor_protocols,2));
    dGMmean_median = zeros(1,size(cor_protocols,2));
    dWMGMdiffmean_boxplot = [];
    dWMGMdiffmean_grp = [];
    dWMGMdiffmean_mean = zeros(1,size(cor_protocols,2));
    dWMGMdiffmean_median = zeros(1,size(cor_protocols,2));
    dWMGMdiffmode_boxplot = [];
    dWMGMdiffmode_grp = [];
    dWMGMdiffmode_mean = zeros(1,size(cor_protocols,2));
    dWMGMdiffmode_median = zeros(1,size(cor_protocols,2));
    dWMmode_boxplot = [];
    dWMmode_grp = [];
    dWMmode_mean = zeros(1,size(cor_protocols,2));
    dWMmode_median = zeros(1,size(cor_protocols,2));
    dGMmode_boxplot = [];
    dGMmode_grp = [];
    dGMmode_mean = zeros(1,size(cor_protocols,2));
    dGMmode_median = zeros(1,size(cor_protocols,2));
    dWMskewness_boxplot = [];
    dWMskewness_grp = [];
    dWMskewness_mean = zeros(1,size(cor_protocols,2));
    dWMskewness_median = zeros(1,size(cor_protocols,2));
    dGMskewness_boxplot = [];
    dGMskewness_grp = [];
    dGMskewness_mean = zeros(1,size(cor_protocols,2));
    dGMskewness_median = zeros(1,size(cor_protocols,2));
    dWMkurtosis_boxplot = [];
    dWMkurtosis_grp = [];
    dWMkurtosis_mean = zeros(1,size(cor_protocols,2));
    dWMkurtosis_median = zeros(1,size(cor_protocols,2));
    dGMkurtosis_boxplot = [];
    dGMkurtosis_grp = [];
    dGMkurtosis_mean = zeros(1,size(cor_protocols,2));
    dGMkurtosis_median = zeros(1,size(cor_protocols,2));
    dwmstd_boxplot = [];
    dwmstd_grp = [];
    dwmstd_mean = zeros(1,size(cor_protocols,2));
    dwmstd_median = zeros(1,size(cor_protocols,2));
    dgmstd_boxplot = [];
    dgmstd_grp = [];
    dgmstd_mean = zeros(1,size(cor_protocols,2));
    dgmstd_median = zeros(1,size(cor_protocols,2));
    fieldcoefmean_boxplot = [];
    fieldcoefmean_grp = [];
    fieldcoefmean_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmean_median = zeros(1,size(cor_protocols,2));
    fieldcoefmedian_boxplot = [];
    fieldcoefmedian_grp = [];
    fieldcoefmedian_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmedian_median = zeros(1,size(cor_protocols,2));
    fieldcoefstd_boxplot = [];
    fieldcoefstd_grp = [];
    fieldcoefstd_mean = zeros(1,size(cor_protocols,2));
    fieldcoefstd_median = zeros(1,size(cor_protocols,2));
    fieldcoefmean3_boxplot = [];
    fieldcoefmean3_grp = [];
    fieldcoefmean3_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmean3_median = zeros(1,size(cor_protocols,2));
    fieldcoefmedian3_boxplot = [];
    fieldcoefmedian3_grp = [];
    fieldcoefmedian3_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmedian3_median = zeros(1,size(cor_protocols,2));
    fieldcoefstd3_boxplot = [];
    fieldcoefstd3_grp = [];
    fieldcoefstd3_mean = zeros(1,size(cor_protocols,2));
    fieldcoefstd3_median = zeros(1,size(cor_protocols,2));
    fieldcoefmean56_boxplot = [];
    fieldcoefmean56_grp = [];
    fieldcoefmean56_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmean56_median = zeros(1,size(cor_protocols,2));
    fieldcoefmedian56_boxplot = [];
    fieldcoefmedian56_grp = [];
    fieldcoefmedian56_mean = zeros(1,size(cor_protocols,2));
    fieldcoefmedian56_median = zeros(1,size(cor_protocols,2));
    fieldcoefstd56_boxplot = [];
    fieldcoefstd56_grp = [];
    fieldcoefstd56_mean = zeros(1,size(cor_protocols,2));
    fieldcoefstd56_median = zeros(1,size(cor_protocols,2));
%% DMRI METRIC POST-PROCESSING INTO FORMAT USABLE FOR SCATTERPLOT VISUALIZATIONS
    for prot = 1:size(cor_protocols,2)
        pom = cor_protocols(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        cor_mean(1,prot)=mean(pom);
        cor_median(1,prot)=median(pom);
        cor_boxplot = [cor_boxplot; pom];
        cor_grp = [cor_grp; pom_grp];

        pom = t_value(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        t_value_mean(1,prot)=mean(pom);
        t_value_median(1,prot)=median(pom);
        t_value_boxplot = [t_value_boxplot; pom];
        t_value_grp = [t_value_grp; pom_grp];

        pom = MI(:,prot);
        pom_grp_sub = [subject{:,2}]';
        pom_grp_sub(pom==0) = [];
        subject_grp = [subject_grp; pom_grp_sub];
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MI_mean(1,prot)=mean(pom);
        MI_median(1,prot)=median(pom);
        MI_boxplot = [MI_boxplot; pom];
        MI_grp = [MI_grp; pom_grp];

        pom = f2_prob(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f2_prob_mean(1,prot)=mean(pom);
        f2_prob_median(1,prot)=median(pom);
        f2_prob_boxplot = [f2_prob_boxplot; pom];
        f2_prob_grp = [f2_prob_grp; pom_grp];

        pom = MIf1(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MIf1_mean(1,prot)=mean(pom);
        MIf1_median(1,prot)=median(pom);
        MIf1_boxplot = [MIf1_boxplot; pom];
        MIf1_grp = [MIf1_grp; pom_grp];

        pom = MIMD(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MIMD_mean(1,prot)=mean(pom);
        MIMD_median(1,prot)=median(pom);
        MIMD_boxplot = [MIMD_boxplot; pom];
        MIMD_grp = [MIMD_grp; pom_grp];

        pom = MId(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MId_mean(1,prot)=mean(pom);
        MId_median(1,prot)=median(pom);
        MId_boxplot = [MId_boxplot; pom];
        MId_grp = [MId_grp; pom_grp];

        pom = FA_wm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWM_mean(1,prot)=mean(pom);
        FAWM_median(1,prot)=median(pom);
        FAWM_boxplot = [FAWM_boxplot; pom];
        FAWM_grp = [FAWM_grp; pom_grp];

        pom = FA_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAGM_mean(1,prot)=mean(pom);
        FAGM_median(1,prot)=median(pom);
        FAGM_boxplot = [FAGM_boxplot; pom];
        FAGM_grp = [FAGM_grp; pom_grp];

        pom = FA_diff_wm_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMGMdiff_mean(1,prot)=mean(pom);
        FAWMGMdiff_median(1,prot)=median(pom);
        FAWMGMdiff_boxplot = [FAWMGMdiff_boxplot; pom];
        FAWMGMdiff_grp = [FAWMGMdiff_grp; pom_grp];

        pom = FA_diff_wm_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMGMdiffmean_mean(1,prot)=mean(pom);
        FAWMGMdiffmean_median(1,prot)=median(pom);
        FAWMGMdiffmean_boxplot = [FAWMGMdiffmean_boxplot; pom];
        FAWMGMdiffmean_grp = [FAWMGMdiffmean_grp; pom_grp];
        
        pom2 = FA_diff_wm_gm_mean(:,prot);
        pom = FA_diff_wm_gm_mode(:,prot);
        pom(pom2==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMGMdiffmode_mean(1,prot)=mean(pom);
        FAWMGMdiffmode_median(1,prot)=median(pom);
        FAWMGMdiffmode_boxplot = [FAWMGMdiffmode_boxplot; pom];
        FAWMGMdiffmode_grp = [FAWMGMdiffmode_grp; pom_grp];

        pom = sqrt(FA_wm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAwmstd_mean(1,prot)=mean(pom);
        FAwmstd_median(1,prot)=median(pom);
        FAwmstd_boxplot = [FAwmstd_boxplot; pom];
        FAwmstd_grp = [FAwmstd_grp; pom_grp];

        pom = sqrt(FA_gm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAgmstd_mean(1,prot)=mean(pom);
        FAgmstd_median(1,prot)=median(pom);
        FAgmstd_boxplot = [FAgmstd_boxplot; pom];
        FAgmstd_grp = [FAgmstd_grp; pom_grp];

        pom = fieldcoef_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmean_mean(1,prot)=mean(pom);
        fieldcoefmean_median(1,prot)=median(pom);
        fieldcoefmean_boxplot = [fieldcoefmean_boxplot; pom];
        fieldcoefmean_grp = [fieldcoefmean_grp; pom_grp];

        pom = fieldcoef_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmedian_mean(1,prot)=mean(pom);
        fieldcoefmedian_median(1,prot)=median(pom);
        fieldcoefmedian_boxplot = [fieldcoefmedian_boxplot; pom];
        fieldcoefmedian_grp = [fieldcoefmedian_grp; pom_grp];

        pom = sqrt(fieldcoef_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefstd_mean(1,prot)=mean(pom);
        fieldcoefstd_median(1,prot)=median(pom);
        fieldcoefstd_boxplot = [fieldcoefstd_boxplot; pom];
        fieldcoefstd_grp = [fieldcoefstd_grp; pom_grp];

        pom = fieldcoef3_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmean3_mean(1,prot)=mean(pom);
        fieldcoefmean3_median(1,prot)=median(pom);
        fieldcoefmean3_boxplot = [fieldcoefmean3_boxplot; pom];
        fieldcoefmean3_grp = [fieldcoefmean3_grp; pom_grp];

        pom = fieldcoef3_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmedian3_mean(1,prot)=mean(pom);
        fieldcoefmedian3_median(1,prot)=median(pom);
        fieldcoefmedian3_boxplot = [fieldcoefmedian3_boxplot; pom];
        fieldcoefmedian3_grp = [fieldcoefmedian3_grp; pom_grp];

        pom = sqrt(fieldcoef3_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefstd3_mean(1,prot)=mean(pom);
        fieldcoefstd3_median(1,prot)=median(pom);
        fieldcoefstd3_boxplot = [fieldcoefstd3_boxplot; pom];
        fieldcoefstd3_grp = [fieldcoefstd3_grp; pom_grp];

        pom = fieldcoef56_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmean56_mean(1,prot)=mean(pom);
        fieldcoefmean56_median(1,prot)=median(pom);
        fieldcoefmean56_boxplot = [fieldcoefmean56_boxplot; pom];
        fieldcoefmean56_grp = [fieldcoefmean56_grp; pom_grp];

        pom = fieldcoef56_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefmedian56_mean(1,prot)=mean(pom);
        fieldcoefmedian56_median(1,prot)=median(pom);
        fieldcoefmedian56_boxplot = [fieldcoefmedian56_boxplot; pom];
        fieldcoefmedian56_grp = [fieldcoefmedian56_grp; pom_grp];

        pom = sqrt(fieldcoef56_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        fieldcoefstd56_mean(1,prot)=mean(pom);
        fieldcoefstd56_median(1,prot)=median(pom);
        fieldcoefstd56_boxplot = [fieldcoefstd56_boxplot; pom];
        fieldcoefstd56_grp = [fieldcoefstd56_grp; pom_grp];

        pom = FA_gm_pKolSmir(:,prot);
        pom(pom==3) = [];
        pom_grp = prot*ones(size(pom));
        FAgmpKolSmir_mean(1,prot)=mean(pom);
        FAgmpKolSmir_median(1,prot)=median(pom);
        FAgmpKolSmir_boxplot = [FAgmpKolSmir_boxplot; pom];
        FAgmpKolSmir_grp = [FAgmpKolSmir_grp; pom_grp];

        pom = FA_wm_pKolSmir(:,prot);
        pom(pom==3) = [];
        pom_grp = prot*ones(size(pom));
        FAwmpKolSmir_mean(1,prot)=mean(pom);
        FAwmpKolSmir_median(1,prot)=median(pom);
        FAwmpKolSmir_boxplot = [FAwmpKolSmir_boxplot; pom];
        FAwmpKolSmir_grp = [FAwmpKolSmir_grp; pom_grp];

        pom = MD_gm_pKolSmir(:,prot);
        pom(pom==3) = [];
        pom_grp = prot*ones(size(pom));
        MDgmpKolSmir_mean(1,prot)=mean(pom);
        MDgmpKolSmir_median(1,prot)=median(pom);
        MDgmpKolSmir_boxplot = [MDgmpKolSmir_boxplot; pom];
        MDgmpKolSmir_grp = [MDgmpKolSmir_grp; pom_grp];

        pom = MD_wm_pKolSmir(:,prot);
        pom(pom==3) = [];
        pom_grp = prot*ones(size(pom));
        MDwmpKolSmir_mean(1,prot)=mean(pom);
        MDwmpKolSmir_median(1,prot)=median(pom);
        MDwmpKolSmir_boxplot = [MDwmpKolSmir_boxplot; pom];
        MDwmpKolSmir_grp = [MDwmpKolSmir_grp; pom_grp];

        pom = MD_gauss_wm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWM_mean(1,prot)=mean(pom);
        MDWM_median(1,prot)=median(pom);
        MDWM_boxplot = [MDWM_boxplot; pom];
        MDWM_grp = [MDWM_grp; pom_grp];

        pom = MD_gauss_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDGM_mean(1,prot)=mean(pom);
        MDGM_median(1,prot)=median(pom);
        MDGM_boxplot = [MDGM_boxplot; pom];
        MDGM_grp = [MDGM_grp; pom_grp];

        pom = MD_gauss_diff_wm_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMGMdiff_mean(1,prot)=mean(pom);
        MDWMGMdiff_median(1,prot)=median(pom);
        MDWMGMdiff_boxplot = [MDWMGMdiff_boxplot; pom];
        MDWMGMdiff_grp = [MDWMGMdiff_grp; pom_grp];

        pom = MD_gauss_diff_wm_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMGMdiffmean_mean(1,prot)=mean(pom);
        MDWMGMdiffmean_median(1,prot)=median(pom);
        MDWMGMdiffmean_boxplot = [MDWMGMdiffmean_boxplot; pom];
        MDWMGMdiffmean_grp = [MDWMGMdiffmean_grp; pom_grp];

        pom2 = MD_gauss_gm_median(:,prot);
        pom = MD_gauss_diff_wm_gm_mode(:,prot);
        pom(pom2==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMGMdiffmode_mean(1,prot)=mean(pom);
        MDWMGMdiffmode_median(1,prot)=median(pom);
        MDWMGMdiffmode_boxplot = [MDWMGMdiffmode_boxplot; pom];
        MDWMGMdiffmode_grp = [MDWMGMdiffmode_grp; pom_grp];

        pom = sqrt(MD_gauss_wm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDwmstd_mean(1,prot)=mean(pom);
        MDwmstd_median(1,prot)=median(pom);
        MDwmstd_boxplot = [MDwmstd_boxplot; pom];
        MDwmstd_grp = [MDwmstd_grp; pom_grp];

        pom = sqrt(MD_gauss_gm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDgmstd_mean(1,prot)=mean(pom);
        MDgmstd_median(1,prot)=median(pom);
        MDgmstd_boxplot = [MDgmstd_boxplot; pom];
        MDgmstd_grp = [MDgmstd_grp; pom_grp];

        pom = f1_gauss_wm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WM_mean(1,prot)=mean(pom);
        f1WM_median(1,prot)=median(pom);
        f1WM_boxplot = [f1WM_boxplot; pom];
        f1WM_grp = [f1WM_grp; pom_grp];

        pom = f1_gauss_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1GM_mean(1,prot)=mean(pom);
        f1GM_median(1,prot)=median(pom);
        f1GM_boxplot = [f1GM_boxplot; pom];
        f1GM_grp = [f1GM_grp; pom_grp];

        pom = f1_gauss_diff_wm_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMGMdiff_mean(1,prot)=mean(pom);
        f1WMGMdiff_median(1,prot)=median(pom);
        f1WMGMdiff_boxplot = [f1WMGMdiff_boxplot; pom];
        f1WMGMdiff_grp = [f1WMGMdiff_grp; pom_grp];

        pom = f1_gauss_diff_wm_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMGMdiffmean_mean(1,prot)=mean(pom);
        f1WMGMdiffmean_median(1,prot)=median(pom);
        f1WMGMdiffmean_boxplot = [f1WMGMdiffmean_boxplot; pom];
        f1WMGMdiffmean_grp = [f1WMGMdiffmean_grp; pom_grp];

        pom2 = f1_gauss_gm_median(:,prot);
        pom = f1_gauss_diff_wm_gm_mode(:,prot);
        pom(pom2==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMGMdiffmode_mean(1,prot)=mean(pom);
        f1WMGMdiffmode_median(1,prot)=median(pom);
        f1WMGMdiffmode_boxplot = [f1WMGMdiffmode_boxplot; pom];
        f1WMGMdiffmode_grp = [f1WMGMdiffmode_grp; pom_grp];

        pom = sqrt(f1_gauss_wm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1wmstd_mean(1,prot)=mean(pom);
        f1wmstd_median(1,prot)=median(pom);
        f1wmstd_boxplot = [f1wmstd_boxplot; pom];
        f1wmstd_grp = [f1wmstd_grp; pom_grp];

        pom = sqrt(f1_gauss_gm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1gmstd_mean(1,prot)=mean(pom);
        f1gmstd_median(1,prot)=median(pom);
        f1gmstd_boxplot = [f1gmstd_boxplot; pom];
        f1gmstd_grp = [f1gmstd_grp; pom_grp];

        pom = d_gauss_wm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWM_mean(1,prot)=mean(pom);
        dWM_median(1,prot)=median(pom);
        dWM_boxplot = [dWM_boxplot; pom];
        dWM_grp = [dWM_grp; pom_grp];

        pom = d_gauss_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dGM_mean(1,prot)=mean(pom);
        dGM_median(1,prot)=median(pom);
        dGM_boxplot = [dGM_boxplot; pom];
        dGM_grp = [dGM_grp; pom_grp];

        pom = d_gauss_diff_wm_gm_median(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMGMdiff_mean(1,prot)=mean(pom);
        dWMGMdiff_median(1,prot)=median(pom);
        dWMGMdiff_boxplot = [dWMGMdiff_boxplot; pom];
        dWMGMdiff_grp = [dWMGMdiff_grp; pom_grp];

        pom = d_gauss_wm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMmean_mean(1,prot)=mean(pom);
        dWMmean_median(1,prot)=median(pom);
        dWMmean_boxplot = [dWMmean_boxplot; pom];
        dWMmean_grp = [dWMmean_grp; pom_grp];

        pom = d_gauss_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dGMmean_mean(1,prot)=mean(pom);
        dGMmean_median(1,prot)=median(pom);
        dGMmean_boxplot = [dGMmean_boxplot; pom];
        dGMmean_grp = [dGMmean_grp; pom_grp];

        pom = d_gauss_diff_wm_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMGMdiffmean_mean(1,prot)=mean(pom);
        dWMGMdiffmean_median(1,prot)=median(pom);
        dWMGMdiffmean_boxplot = [dWMGMdiffmean_boxplot; pom];
        dWMGMdiffmean_grp = [dWMGMdiffmean_grp; pom_grp];

        pom2 = d_gauss_wm_mean(:,prot);
        pom = d_gauss_diff_wm_gm_mode(:,prot);
        pom(pom2==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMGMdiffmode_mean(1,prot)=mean(pom);
        dWMGMdiffmode_median(1,prot)=median(pom);
        dWMGMdiffmode_boxplot = [dWMGMdiffmode_boxplot; pom];
        dWMGMdiffmode_grp = [dWMGMdiffmode_grp; pom_grp];

        pom = d_gauss_wm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMmode_mean(1,prot)=mean(pom);
        dWMmode_median(1,prot)=median(pom);
        dWMmode_boxplot = [dWMmode_boxplot; pom];
        dWMmode_grp = [dWMmode_grp; pom_grp];

        pom = d_gauss_gm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dGMmode_mean(1,prot)=mean(pom);
        dGMmode_median(1,prot)=median(pom);
        dGMmode_boxplot = [dGMmode_boxplot; pom];
        dGMmode_grp = [dGMmode_grp; pom_grp];

        pom = sqrt(d_gauss_wm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dwmstd_mean(1,prot)=mean(pom);
        dwmstd_median(1,prot)=median(pom);
        dwmstd_boxplot = [dwmstd_boxplot; pom];
        dwmstd_grp = [dwmstd_grp; pom_grp];

        pom = sqrt(d_gauss_gm_var(:,prot));
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dgmstd_mean(1,prot)=mean(pom);
        dgmstd_median(1,prot)=median(pom);
        dgmstd_boxplot = [dgmstd_boxplot; pom];
        dgmstd_grp = [dgmstd_grp; pom_grp];

        pom = FA_wm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMmean_mean(1,prot)=mean(pom);
        FAWMmean_median(1,prot)=median(pom);
        FAWMmean_boxplot = [FAWMmean_boxplot; pom];
        FAWMmean_grp = [FAWMmean_grp; pom_grp];

        pom = FA_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAGMmean_mean(1,prot)=mean(pom);
        FAGMmean_median(1,prot)=median(pom);
        FAGMmean_boxplot = [FAGMmean_boxplot; pom];
        FAGMmean_grp = [FAGMmean_grp; pom_grp];

        pom = FA_wm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMmode_mean(1,prot)=mean(pom);
        FAWMmode_median(1,prot)=median(pom);
        FAWMmode_boxplot = [FAWMmode_boxplot; pom];
        FAWMmode_grp = [FAWMmode_grp; pom_grp];

        pom = FA_gm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAGMmode_mean(1,prot)=mean(pom);
        FAGMmode_median(1,prot)=median(pom);
        FAGMmode_boxplot = [FAGMmode_boxplot; pom];
        FAGMmode_grp = [FAGMmode_grp; pom_grp];

        pom = f1_gauss_wm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMmean_mean(1,prot)=mean(pom);
        f1WMmean_median(1,prot)=median(pom);
        f1WMmean_boxplot = [f1WMmean_boxplot; pom];
        f1WMmean_grp = [f1WMmean_grp; pom_grp];

        pom = f1_gauss_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1GMmean_mean(1,prot)=mean(pom);
        f1GMmean_median(1,prot)=median(pom);
        f1GMmean_boxplot = [f1GMmean_boxplot; pom];
        f1GMmean_grp = [f1GMmean_grp; pom_grp];

        pom = f1_gauss_wm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMmode_mean(1,prot)=mean(pom);
        f1WMmode_median(1,prot)=median(pom);
        f1WMmode_boxplot = [f1WMmode_boxplot; pom];
        f1WMmode_grp = [f1WMmode_grp; pom_grp];

        pom = f1_gauss_gm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1GMmode_mean(1,prot)=mean(pom);
        f1GMmode_median(1,prot)=median(pom);
        f1GMmode_boxplot = [f1GMmode_boxplot; pom];
        f1GMmode_grp = [f1GMmode_grp; pom_grp];

        pom = MD_gauss_wm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMmean_mean(1,prot)=mean(pom);
        MDWMmean_median(1,prot)=median(pom);
        MDWMmean_boxplot = [MDWMmean_boxplot; pom];
        MDWMmean_grp = [MDWMmean_grp; pom_grp];

        pom = MD_gauss_gm_mean(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDGMmean_mean(1,prot)=mean(pom);
        MDGMmean_median(1,prot)=median(pom);
        MDGMmean_boxplot = [MDGMmean_boxplot; pom];
        MDGMmean_grp = [MDGMmean_grp; pom_grp];

        pom = MD_gauss_wm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMmode_mean(1,prot)=mean(pom);
        MDWMmode_median(1,prot)=median(pom);
        MDWMmode_boxplot = [MDWMmode_boxplot; pom];
        MDWMmode_grp = [MDWMmode_grp; pom_grp];

        pom = MD_gauss_gm_mode(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDGMmode_mean(1,prot)=mean(pom);
        MDGMmode_median(1,prot)=median(pom);
        MDGMmode_boxplot = [MDGMmode_boxplot; pom];
        MDGMmode_grp = [MDGMmode_grp; pom_grp];

        pom = FA_gm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAGMskewness_mean(1,prot)=mean(pom);
        FAGMskewness_median(1,prot)=median(pom);
        FAGMskewness_boxplot = [FAGMskewness_boxplot; pom];
        FAGMskewness_grp = [FAGMskewness_grp; pom_grp];

        pom = FA_wm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMskewness_mean(1,prot)=mean(pom);
        FAWMskewness_median(1,prot)=median(pom);
        FAWMskewness_boxplot = [FAWMskewness_boxplot; pom];
        FAWMskewness_grp = [FAWMskewness_grp; pom_grp];

        pom = FA_gm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAGMkurtosis_mean(1,prot)=mean(pom);
        FAGMkurtosis_median(1,prot)=median(pom);
        FAGMkurtosis_boxplot = [FAGMkurtosis_boxplot; pom];
        FAGMkurtosis_grp = [FAGMkurtosis_grp; pom_grp];

        pom = FA_wm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        FAWMkurtosis_mean(1,prot)=mean(pom);
        FAWMkurtosis_median(1,prot)=median(pom);
        FAWMkurtosis_boxplot = [FAWMkurtosis_boxplot; pom];
        FAWMkurtosis_grp = [FAWMkurtosis_grp; pom_grp];

        pom = f1_gauss_gm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1GMskewness_mean(1,prot)=mean(pom);
        f1GMskewness_median(1,prot)=median(pom);
        f1GMskewness_boxplot = [f1GMskewness_boxplot; pom];
        f1GMskewness_grp = [f1GMskewness_grp; pom_grp];

        pom = f1_gauss_wm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMskewness_mean(1,prot)=mean(pom);
        f1WMskewness_median(1,prot)=median(pom);
        f1WMskewness_boxplot = [f1WMskewness_boxplot; pom];
        f1WMskewness_grp = [f1WMskewness_grp; pom_grp];

        pom = f1_gauss_gm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1GMkurtosis_mean(1,prot)=mean(pom);
        f1GMkurtosis_median(1,prot)=median(pom);
        f1GMkurtosis_boxplot = [f1GMkurtosis_boxplot; pom];
        f1GMkurtosis_grp = [f1GMkurtosis_grp; pom_grp];

        pom = f1_gauss_wm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        f1WMkurtosis_mean(1,prot)=mean(pom);
        f1WMkurtosis_median(1,prot)=median(pom);
        f1WMkurtosis_boxplot = [f1WMkurtosis_boxplot; pom];
        f1WMkurtosis_grp = [f1WMkurtosis_grp; pom_grp];

        pom = MD_gauss_gm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDGMskewness_mean(1,prot)=mean(pom);
        MDGMskewness_median(1,prot)=median(pom);
        MDGMskewness_boxplot = [MDGMskewness_boxplot; pom];
        MDGMskewness_grp = [MDGMskewness_grp; pom_grp];

        pom = MD_gauss_wm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMskewness_mean(1,prot)=mean(pom);
        MDWMskewness_median(1,prot)=median(pom);
        MDWMskewness_boxplot = [MDWMskewness_boxplot; pom];
        MDWMskewness_grp = [MDWMskewness_grp; pom_grp];

        pom = MD_gauss_gm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDGMkurtosis_mean(1,prot)=mean(pom);
        MDGMkurtosis_median(1,prot)=median(pom);
        MDGMkurtosis_boxplot = [MDGMkurtosis_boxplot; pom];
        MDGMkurtosis_grp = [MDGMkurtosis_grp; pom_grp];

        pom = MD_gauss_wm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        MDWMkurtosis_mean(1,prot)=mean(pom);
        MDWMkurtosis_median(1,prot)=median(pom);
        MDWMkurtosis_boxplot = [MDWMkurtosis_boxplot; pom];
        MDWMkurtosis_grp = [MDWMkurtosis_grp; pom_grp];

        pom = d_gauss_gm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dGMskewness_mean(1,prot)=mean(pom);
        dGMskewness_median(1,prot)=median(pom);
        dGMskewness_boxplot = [dGMskewness_boxplot; pom];
        dGMskewness_grp = [dGMskewness_grp; pom_grp];

        pom = d_gauss_wm_skewness(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMskewness_mean(1,prot)=mean(pom);
        dWMskewness_median(1,prot)=median(pom);
        dWMskewness_boxplot = [dWMskewness_boxplot; pom];
        dWMskewness_grp = [dWMskewness_grp; pom_grp];

        pom = d_gauss_gm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dGMkurtosis_mean(1,prot)=mean(pom);
        dGMkurtosis_median(1,prot)=median(pom);
        dGMkurtosis_boxplot = [dGMkurtosis_boxplot; pom];
        dGMkurtosis_grp = [dGMkurtosis_grp; pom_grp];

        pom = d_gauss_wm_kurtosis(:,prot);
        pom(pom==0) = [];
        pom_grp = prot*ones(size(pom));
        dWMkurtosis_mean(1,prot)=mean(pom);
        dWMkurtosis_median(1,prot)=median(pom);
        dWMkurtosis_boxplot = [dWMkurtosis_boxplot; pom];
        dWMkurtosis_grp = [dWMkurtosis_grp; pom_grp];

        disp(['Post-processing session ' num2str(prot) ' done.'])
    end
    save(workspace_file)
else
    plot_results2 = plot_results;
    load(workspace_file)
    plot_results = plot_results2;
end

if plot_results == 1
    %% VISUALIZATION OF NORMALIZED MUTUTAL INFORMATION BETWEEN FA MAP AND T2* IMAGE WITH WM/GM CONTRAST 
    pCR_ZOOMint = ranksum(MI_boxplot(subject_grp==1 & MI_grp==1),MI_boxplot(subject_grp==3 & MI_grp==1));
    pCP_ZOOMint = ranksum(MI_boxplot(subject_grp==1 & MI_grp==1),MI_boxplot(subject_grp==2 & MI_grp==1));
    pRP_ZOOMint = ranksum(MI_boxplot(subject_grp==3 & MI_grp==1),MI_boxplot(subject_grp==2 & MI_grp==1));
    pSC_ZOOMint = ranksum(MI_boxplot(subject_grp==4 & MI_grp==1),MI_boxplot(subject_grp==1 & MI_grp==1));
    pSR_ZOOMint = ranksum(MI_boxplot(subject_grp==3 & MI_grp==1),MI_boxplot(subject_grp==4 & MI_grp==1));
    pSP_ZOOMint = ranksum(MI_boxplot(subject_grp==4 & MI_grp==1),MI_boxplot(subject_grp==2 & MI_grp==1));
    pCR_ZOOMnotint = ranksum(MI_boxplot(subject_grp==1 & MI_grp==2),MI_boxplot(subject_grp==3 & MI_grp==2));
    pCP_ZOOMnotint = ranksum(MI_boxplot(subject_grp==1 & MI_grp==2),MI_boxplot(subject_grp==2 & MI_grp==2));
    pRP_ZOOMnotint = ranksum(MI_boxplot(subject_grp==3 & MI_grp==2),MI_boxplot(subject_grp==2 & MI_grp==2));
    pSC_ZOOMnotint = ranksum(MI_boxplot(subject_grp==4 & MI_grp==2),MI_boxplot(subject_grp==1 & MI_grp==2));
    pSR_ZOOMnotint = ranksum(MI_boxplot(subject_grp==3 & MI_grp==2),MI_boxplot(subject_grp==4 & MI_grp==2));
    pSP_ZOOMnotint = ranksum(MI_boxplot(subject_grp==4 & MI_grp==2),MI_boxplot(subject_grp==2 & MI_grp==2));
    pCR_RESOLVE = ranksum(MI_boxplot(subject_grp==1 & MI_grp==3),MI_boxplot(subject_grp==3 & MI_grp==3));
    pCP_RESOLVE = ranksum(MI_boxplot(subject_grp==1 & MI_grp==3),MI_boxplot(subject_grp==2 & MI_grp==3));
    pRP_RESOLVE = ranksum(MI_boxplot(subject_grp==3 & MI_grp==3),MI_boxplot(subject_grp==2 & MI_grp==3));
    pSC_RESOLVE = ranksum(MI_boxplot(subject_grp==4 & MI_grp==3),MI_boxplot(subject_grp==1 & MI_grp==3));
    pSR_RESOLVE = ranksum(MI_boxplot(subject_grp==3 & MI_grp==3),MI_boxplot(subject_grp==4 & MI_grp==3));
    pSP_RESOLVE = ranksum(MI_boxplot(subject_grp==4 & MI_grp==3),MI_boxplot(subject_grp==2 & MI_grp==3));
    % p-values of between-group differences are stored in table_pvals variable
    table_pvals(1,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(2).fig = figure(2);
    set(h(2).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MI_grp)/2-0.20 unique(MI_grp)/2+0.20]', repmat(MI_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MI_grp)/2-0.20 unique(MI_grp)/2+0.20]', repmat(MI_median',1,2)','m-','LineWidth',8)
    scatter(MI_grp(subject_grp==3)/2, MI_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MI_grp(subject_grp==1)/2, MI_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MI_grp(subject_grp==4)/2, MI_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MI_grp(subject_grp==2)/2, MI_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.967,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.963,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.959,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.955,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.951,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.947,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Similarity measurements in T_2 axial space between';'FA maps and T_2 axial scans from C3-C6 segments'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Normalized mutual information','FontSize',18)
    axis([0.5/2 3.5/2 0.90 0.97])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mutual_information_over_protocols']), '-dpng', '-r300')

    %% VISUALIZATION OF f2 PERCENTILES 
    pCR_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==3 & f2_prob_grp==1));
    pCP_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==1));
    pRP_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==1));
    pSC_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==1 & f2_prob_grp==1));
    pSR_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==4 & f2_prob_grp==1));
    pSP_ZOOMint = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==1),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==1));
    pCR_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==3 & f2_prob_grp==2));
    pCP_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==2));
    pRP_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==2));
    pSC_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==1 & f2_prob_grp==2));
    pSR_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==4 & f2_prob_grp==2));
    pSP_ZOOMnotint = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==2),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==2));
    pCR_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==3 & f2_prob_grp==3));
    pCP_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==1 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==3));
    pRP_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==3));
    pSC_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==1 & f2_prob_grp==3));
    pSR_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==3 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==4 & f2_prob_grp==3));
    pSP_RESOLVE = ranksum(f2_prob_boxplot(subject_grp==4 & f2_prob_grp==3),f2_prob_boxplot(subject_grp==2 & f2_prob_grp==3));
    % p-values of between-group differences are stored in table_pvals variable
    table_pvals(2,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(3).fig = figure(3);
    set(h(3).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f2_prob_grp)/2-0.20 unique(f2_prob_grp)/2+0.20]', repmat(f2_prob_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f2_prob_grp)/2-0.20 unique(f2_prob_grp)/2+0.20]', repmat(f2_prob_median',1,2)','m-','LineWidth',8)
    scatter(f2_prob_grp(subject_grp==3)/2, f2_prob_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f2_prob_grp(subject_grp==1)/2, f2_prob_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f2_prob_grp(subject_grp==4)/2, f2_prob_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f2_prob_grp(subject_grp==2)/2, f2_prob_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,22.5,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,22.5,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,21.8,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,21.8,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,21.1,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,21.1,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,20.4,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,20.4,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,19.7,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,19.7,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,19.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,19.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,22.5,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,22.5,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,21.8,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,21.8,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,21.1,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,21.1,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,20.4,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,20.4,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,19.7,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,19.7,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,19.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,19.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,22.5,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,22.5,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,21.8,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,21.8,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,21.1,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,21.1,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,20.4,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,20.4,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,19.7,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,19.7,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,19.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,19.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Percentile of spinal cord voxels with volume f_2 > 0.05';'estimated from C3-C6 segments in T_2 axial space'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Percentile [%]','FontSize',18)
    axis([0.5/2 3.5/2 floor(0.98*min(f2_prob_boxplot)) ceil(1.15*max(f2_prob_boxplot)) ])
    % axis([0.25 1.75 0.90 0.97])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f2_percentile_2nd_direction']), '-dpng', '-r300')
    
    %% VISUALIZATION OF NORMALIZED MUTUTAL INFORMATION BETWEEN f1 MAP AND T2* IMAGE WITH WM/GM CONTRAST 
    pCR_ZOOMint = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==1),MIf1_boxplot(subject_grp==3 & MIf1_grp==1));
    pCP_ZOOMint = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==1),MIf1_boxplot(subject_grp==2 & MIf1_grp==1));
    pRP_ZOOMint = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==1),MIf1_boxplot(subject_grp==2 & MIf1_grp==1));
    pSC_ZOOMint = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==1),MIf1_boxplot(subject_grp==1 & MIf1_grp==1));
    pSR_ZOOMint = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==1),MIf1_boxplot(subject_grp==4 & MIf1_grp==1));
    pSP_ZOOMint = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==1),MIf1_boxplot(subject_grp==2 & MIf1_grp==1));
    pCR_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==2),MIf1_boxplot(subject_grp==3 & MIf1_grp==2));
    pCP_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==2),MIf1_boxplot(subject_grp==2 & MIf1_grp==2));
    pRP_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==2),MIf1_boxplot(subject_grp==2 & MIf1_grp==2));
    pSC_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==2),MIf1_boxplot(subject_grp==1 & MIf1_grp==2));
    pSR_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==2),MIf1_boxplot(subject_grp==4 & MIf1_grp==2));
    pSP_ZOOMnotint = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==2),MIf1_boxplot(subject_grp==2 & MIf1_grp==2));
    pCR_RESOLVE = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==3),MIf1_boxplot(subject_grp==3 & MIf1_grp==3));
    pCP_RESOLVE = ranksum(MIf1_boxplot(subject_grp==1 & MIf1_grp==3),MIf1_boxplot(subject_grp==2 & MIf1_grp==3));
    pRP_RESOLVE = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==3),MIf1_boxplot(subject_grp==2 & MIf1_grp==3));
    pSC_RESOLVE = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==3),MIf1_boxplot(subject_grp==1 & MIf1_grp==3));
    pSR_RESOLVE = ranksum(MIf1_boxplot(subject_grp==3 & MIf1_grp==3),MIf1_boxplot(subject_grp==4 & MIf1_grp==3));
    pSP_RESOLVE = ranksum(MIf1_boxplot(subject_grp==4 & MIf1_grp==3),MIf1_boxplot(subject_grp==2 & MIf1_grp==3));
    % p-values of between-group differences are stored in table_pvals variable
    table_pvals(3,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(4).fig = figure(4);
    set(h(4).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MIf1_grp)/2-0.20 unique(MIf1_grp)/2+0.20]', repmat(MIf1_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MIf1_grp)/2-0.20 unique(MIf1_grp)/2+0.20]', repmat(MIf1_median',1,2)','m-','LineWidth',8)
    scatter(MIf1_grp(subject_grp==3)/2, MIf1_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MIf1_grp(subject_grp==1)/2, MIf1_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MIf1_grp(subject_grp==4)/2, MIf1_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MIf1_grp(subject_grp==2)/2, MIf1_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.967,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.963,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.959,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.955,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.951,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.947,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Similarity measurements in T_2 axial space between';'f_1 maps and T_2 axial scans from C3-C6 segments'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Normalized mutual information','FontSize',18)
    axis([0.5/2 3.5/2 0.90 0.97])
    % axis([0.25 1.75 0.90 0.97])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mutual_information_over_protocols']), '-dpng', '-r300')
    
    %% VISUALIZATION OF NORMALIZED MUTUTAL INFORMATION BETWEEN MD MAP AND T2* IMAGE WITH WM/GM CONTRAST 
    pCR_ZOOMint = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==1),MIMD_boxplot(subject_grp==3 & MIMD_grp==1));
    pCP_ZOOMint = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==1),MIMD_boxplot(subject_grp==2 & MIMD_grp==1));
    pRP_ZOOMint = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==1),MIMD_boxplot(subject_grp==2 & MIMD_grp==1));
    pSC_ZOOMint = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==1),MIMD_boxplot(subject_grp==1 & MIMD_grp==1));
    pSR_ZOOMint = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==1),MIMD_boxplot(subject_grp==4 & MIMD_grp==1));
    pSP_ZOOMint = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==1),MIMD_boxplot(subject_grp==2 & MIMD_grp==1));
    pCR_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==2),MIMD_boxplot(subject_grp==3 & MIMD_grp==2));
    pCP_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==2),MIMD_boxplot(subject_grp==2 & MIMD_grp==2));
    pRP_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==2),MIMD_boxplot(subject_grp==2 & MIMD_grp==2));
    pSC_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==2),MIMD_boxplot(subject_grp==1 & MIMD_grp==2));
    pSR_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==2),MIMD_boxplot(subject_grp==4 & MIMD_grp==2));
    pSP_ZOOMnotint = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==2),MIMD_boxplot(subject_grp==2 & MIMD_grp==2));
    pCR_RESOLVE = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==3),MIMD_boxplot(subject_grp==3 & MIMD_grp==3));
    pCP_RESOLVE = ranksum(MIMD_boxplot(subject_grp==1 & MIMD_grp==3),MIMD_boxplot(subject_grp==2 & MIMD_grp==3));
    pRP_RESOLVE = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==3),MIMD_boxplot(subject_grp==2 & MIMD_grp==3));
    pSC_RESOLVE = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==3),MIMD_boxplot(subject_grp==1 & MIMD_grp==3));
    pSR_RESOLVE = ranksum(MIMD_boxplot(subject_grp==3 & MIMD_grp==3),MIMD_boxplot(subject_grp==4 & MIMD_grp==3));
    pSP_RESOLVE = ranksum(MIMD_boxplot(subject_grp==4 & MIMD_grp==3),MIMD_boxplot(subject_grp==2 & MIMD_grp==3));
    table_pvals(4,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(5).fig = figure(5);
    set(h(5).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MIMD_grp)/2-0.20 unique(MIMD_grp)/2+0.20]', repmat(MIMD_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MIMD_grp)/2-0.20 unique(MIMD_grp)/2+0.20]', repmat(MIMD_median',1,2)','m-','LineWidth',8)
    scatter(MIMD_grp(subject_grp==3)/2, MIMD_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MIMD_grp(subject_grp==1)/2, MIMD_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MIMD_grp(subject_grp==4)/2, MIMD_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MIMD_grp(subject_grp==2)/2, MIMD_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.967,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.963,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.959,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.955,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.951,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.947,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Similarity measurements in T_2 axial space between';'MD maps and T_2 axial scans from C3-C6 segments'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Normalized mutual information','FontSize',18)
    axis([0.5/2 3.5/2 0.90 0.97])
    % axis([0.25 1.75 0.90 0.97])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mutual_information_over_protocols']), '-dpng', '-r300')

    %% VISUALIZATION OF NORMALIZED MUTUTAL INFORMATION BETWEEN d MAP AND T2* IMAGE WITH WM/GM CONTRAST 
    pCR_ZOOMint = ranksum(MId_boxplot(subject_grp==1 & MId_grp==1),MId_boxplot(subject_grp==3 & MId_grp==1));
    pCP_ZOOMint = ranksum(MId_boxplot(subject_grp==1 & MId_grp==1),MId_boxplot(subject_grp==2 & MId_grp==1));
    pRP_ZOOMint = ranksum(MId_boxplot(subject_grp==3 & MId_grp==1),MId_boxplot(subject_grp==2 & MId_grp==1));
    pSC_ZOOMint = ranksum(MId_boxplot(subject_grp==4 & MId_grp==1),MId_boxplot(subject_grp==1 & MId_grp==1));
    pSR_ZOOMint = ranksum(MId_boxplot(subject_grp==3 & MId_grp==1),MId_boxplot(subject_grp==4 & MId_grp==1));
    pSP_ZOOMint = ranksum(MId_boxplot(subject_grp==4 & MId_grp==1),MId_boxplot(subject_grp==2 & MId_grp==1));
    pCR_ZOOMnotint = ranksum(MId_boxplot(subject_grp==1 & MId_grp==2),MId_boxplot(subject_grp==3 & MId_grp==2));
    pCP_ZOOMnotint = ranksum(MId_boxplot(subject_grp==1 & MId_grp==2),MId_boxplot(subject_grp==2 & MId_grp==2));
    pRP_ZOOMnotint = ranksum(MId_boxplot(subject_grp==3 & MId_grp==2),MId_boxplot(subject_grp==2 & MId_grp==2));
    pSC_ZOOMnotint = ranksum(MId_boxplot(subject_grp==4 & MId_grp==2),MId_boxplot(subject_grp==1 & MId_grp==2));
    pSR_ZOOMnotint = ranksum(MId_boxplot(subject_grp==3 & MId_grp==2),MId_boxplot(subject_grp==4 & MId_grp==2));
    pSP_ZOOMnotint = ranksum(MId_boxplot(subject_grp==4 & MId_grp==2),MId_boxplot(subject_grp==2 & MId_grp==2));
    pCR_RESOLVE = ranksum(MId_boxplot(subject_grp==1 & MId_grp==3),MId_boxplot(subject_grp==3 & MId_grp==3));
    pCP_RESOLVE = ranksum(MId_boxplot(subject_grp==1 & MId_grp==3),MId_boxplot(subject_grp==2 & MId_grp==3));
    pRP_RESOLVE = ranksum(MId_boxplot(subject_grp==3 & MId_grp==3),MId_boxplot(subject_grp==2 & MId_grp==3));
    pSC_RESOLVE = ranksum(MId_boxplot(subject_grp==4 & MId_grp==3),MId_boxplot(subject_grp==1 & MId_grp==3));
    pSR_RESOLVE = ranksum(MId_boxplot(subject_grp==3 & MId_grp==3),MId_boxplot(subject_grp==4 & MId_grp==3));
    pSP_RESOLVE = ranksum(MId_boxplot(subject_grp==4 & MId_grp==3),MId_boxplot(subject_grp==2 & MId_grp==3));
    table_pvals(5,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(6).fig = figure(6);
    set(h(6).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MId_grp)/2-0.20 unique(MId_grp)/2+0.20]', repmat(MId_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MId_grp)/2-0.20 unique(MId_grp)/2+0.20]', repmat(MId_median',1,2)','m-','LineWidth',8)
    scatter(MId_grp(subject_grp==3)/2, MId_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MId_grp(subject_grp==1)/2, MId_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MId_grp(subject_grp==4)/2, MId_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MId_grp(subject_grp==2)/2, MId_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.967,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.963,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.959,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.955,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.951,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.947,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.967,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.967,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.963,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.963,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.959,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.959,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.955,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.955,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.951,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.951,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.947,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.947,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Similarity measurements in T_2 axial space between mean';'diffusivity maps and T_2 axial scans from C3-C6 segments'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Normalized mutual information','FontSize',18)
    axis([0.5/2 3.5/2 0.90 0.97])
    % axis([0.25 1.75 0.90 0.97])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_meanD_mutual_information_over_protocols']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEDIAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==1),FAWM_boxplot(subject_grp==3 & FAWM_grp==1));
    pCP_ZOOMint = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==1),FAWM_boxplot(subject_grp==2 & FAWM_grp==1));
    pRP_ZOOMint = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==1),FAWM_boxplot(subject_grp==2 & FAWM_grp==1));
    pSC_ZOOMint = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==1),FAWM_boxplot(subject_grp==1 & FAWM_grp==1));
    pSR_ZOOMint = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==1),FAWM_boxplot(subject_grp==4 & FAWM_grp==1));
    pSP_ZOOMint = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==1),FAWM_boxplot(subject_grp==2 & FAWM_grp==1));
    pCR_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==2),FAWM_boxplot(subject_grp==3 & FAWM_grp==2));
    pCP_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==2),FAWM_boxplot(subject_grp==2 & FAWM_grp==2));
    pRP_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==2),FAWM_boxplot(subject_grp==2 & FAWM_grp==2));
    pSC_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==2),FAWM_boxplot(subject_grp==1 & FAWM_grp==2));
    pSR_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==2),FAWM_boxplot(subject_grp==4 & FAWM_grp==2));
    pSP_ZOOMnotint = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==2),FAWM_boxplot(subject_grp==2 & FAWM_grp==2));
    pCR_RESOLVE = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==3),FAWM_boxplot(subject_grp==3 & FAWM_grp==3));
    pCP_RESOLVE = ranksum(FAWM_boxplot(subject_grp==1 & FAWM_grp==3),FAWM_boxplot(subject_grp==2 & FAWM_grp==3));
    pRP_RESOLVE = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==3),FAWM_boxplot(subject_grp==2 & FAWM_grp==3));
    pSC_RESOLVE = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==3),FAWM_boxplot(subject_grp==1 & FAWM_grp==3));
    pSR_RESOLVE = ranksum(FAWM_boxplot(subject_grp==3 & FAWM_grp==3),FAWM_boxplot(subject_grp==4 & FAWM_grp==3));
    pSP_RESOLVE = ranksum(FAWM_boxplot(subject_grp==4 & FAWM_grp==3),FAWM_boxplot(subject_grp==2 & FAWM_grp==3));
    table_pvals(6,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(7).fig = figure(7);
    set(h(7).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWM_grp)/2-0.20 unique(FAWM_grp)/2+0.20]', repmat(FAWM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWM_grp)/2-0.20 unique(FAWM_grp)/2+0.20]', repmat(FAWM_median',1,2)','m-','LineWidth',8)
    scatter(FAWM_grp(subject_grp==3)/2, FAWM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWM_grp(subject_grp==1)/2, FAWM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWM_grp(subject_grp==4)/2, FAWM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWM_grp(subject_grp==2)/2, FAWM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.739,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.721,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.712,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.703,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.694,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA median value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA median values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.80])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_median_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEDIAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==1),FAGM_boxplot(subject_grp==3 & FAGM_grp==1));
    pCP_ZOOMint = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==1),FAGM_boxplot(subject_grp==2 & FAGM_grp==1));
    pRP_ZOOMint = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==1),FAGM_boxplot(subject_grp==2 & FAGM_grp==1));
    pSC_ZOOMint = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==1),FAGM_boxplot(subject_grp==1 & FAGM_grp==1));
    pSR_ZOOMint = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==1),FAGM_boxplot(subject_grp==4 & FAGM_grp==1));
    pSP_ZOOMint = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==1),FAGM_boxplot(subject_grp==2 & FAGM_grp==1));
    pCR_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==2),FAGM_boxplot(subject_grp==3 & FAGM_grp==2));
    pCP_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==2),FAGM_boxplot(subject_grp==2 & FAGM_grp==2));
    pRP_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==2),FAGM_boxplot(subject_grp==2 & FAGM_grp==2));
    pSC_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==2),FAGM_boxplot(subject_grp==1 & FAGM_grp==2));
    pSR_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==2),FAGM_boxplot(subject_grp==4 & FAGM_grp==2));
    pSP_ZOOMnotint = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==2),FAGM_boxplot(subject_grp==2 & FAGM_grp==2));
    pCR_RESOLVE = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==3),FAGM_boxplot(subject_grp==3 & FAGM_grp==3));
    pCP_RESOLVE = ranksum(FAGM_boxplot(subject_grp==1 & FAGM_grp==3),FAGM_boxplot(subject_grp==2 & FAGM_grp==3));
    pRP_RESOLVE = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==3),FAGM_boxplot(subject_grp==2 & FAGM_grp==3));
    pSC_RESOLVE = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==3),FAGM_boxplot(subject_grp==1 & FAGM_grp==3));
    pSR_RESOLVE = ranksum(FAGM_boxplot(subject_grp==3 & FAGM_grp==3),FAGM_boxplot(subject_grp==4 & FAGM_grp==3));
    pSP_RESOLVE = ranksum(FAGM_boxplot(subject_grp==4 & FAGM_grp==3),FAGM_boxplot(subject_grp==2 & FAGM_grp==3));
    table_pvals(7,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(8).fig = figure(8);
    set(h(8).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAGM_grp)/2-0.20 unique(FAGM_grp)/2+0.20]', repmat(FAGM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAGM_grp)/2-0.20 unique(FAGM_grp)/2+0.20]', repmat(FAGM_median',1,2)','m-','LineWidth',8)
    scatter(FAGM_grp(subject_grp==3)/2, FAGM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGM_grp(subject_grp==1)/2, FAGM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGM_grp(subject_grp==4)/2, FAGM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAGM_grp(subject_grp==2)/2, FAGM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.739,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.721,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.712,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.703,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.694,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA median value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA median values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.80])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_median_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEDIAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==1));
    pCP_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==1));
    pRP_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==1));
    pSC_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==1));
    pSR_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==1));
    pSP_ZOOMint = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==1),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==2),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==2));
    pCR_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==3));
    pCP_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==3));
    pRP_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==3));
    pSC_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==1 & FAWMGMdiff_grp==3));
    pSR_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==3 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==3));
    pSP_RESOLVE = ranksum(FAWMGMdiff_boxplot(subject_grp==4 & FAWMGMdiff_grp==3),FAWMGMdiff_boxplot(subject_grp==2 & FAWMGMdiff_grp==3));
    table_pvals(8,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(9).fig = figure(9);
    set(h(9).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMGMdiff_grp)/2-0.20 unique(FAWMGMdiff_grp)/2+0.20]', repmat(FAWMGMdiff_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMGMdiff_grp)/2-0.20 unique(FAWMGMdiff_grp)/2+0.20]', repmat(FAWMGMdiff_median',1,2)','m-','LineWidth',8)
    scatter(FAWMGMdiff_grp(subject_grp==3)/2, FAWMGMdiff_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiff_grp(subject_grp==1)/2, FAWMGMdiff_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiff_grp(subject_grp==4)/2, FAWMGMdiff_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiff_grp(subject_grp==2)/2, FAWMGMdiff_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.18,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.16,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.14,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.12,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.08,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'median FA values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA median differences','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_median_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA STD VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==1));
    pCP_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==1));
    pRP_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==1));
    pSC_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==1));
    pSR_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==1));
    pSP_ZOOMint = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==1),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==1));
    pCR_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==2));
    pCP_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==2));
    pRP_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==2));
    pSC_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==2));
    pSR_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==2));
    pSP_ZOOMnotint = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==2),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==2));
    pCR_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==3));
    pCP_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==3));
    pRP_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==3));
    pSC_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==1 & FAwmstd_grp==3));
    pSR_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==3 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==3));
    pSP_RESOLVE = ranksum(FAwmstd_boxplot(subject_grp==4 & FAwmstd_grp==3),FAwmstd_boxplot(subject_grp==2 & FAwmstd_grp==3));
    table_pvals(9,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(10).fig = figure(10);
    set(h(10).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAwmstd_grp)/2-0.20 unique(FAwmstd_grp)/2+0.20]', repmat(FAwmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAwmstd_grp)/2-0.20 unique(FAwmstd_grp)/2+0.20]', repmat(FAwmstd_median',1,2)','m-','LineWidth',8)
    scatter(FAwmstd_grp(subject_grp==3)/2, FAwmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAwmstd_grp(subject_grp==1)/2, FAwmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAwmstd_grp(subject_grp==4)/2, FAwmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAwmstd_grp(subject_grp==2)/2, FAwmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.185,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.180,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.175,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.170,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.165,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.160,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.185,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.180,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.175,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.170,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.165,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.160,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.185,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.180,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.175,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.170,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.165,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.160,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA standard deviation value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.06 0.190])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_std_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA STD VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==1));
    pCP_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==1));
    pRP_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==1));
    pSC_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==1));
    pSR_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==1));
    pSP_ZOOMint = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==1),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==1));
    pCR_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==2));
    pCP_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==2));
    pRP_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==2));
    pSC_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==2));
    pSR_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==2));
    pSP_ZOOMnotint = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==2),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==2));
    pCR_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==3));
    pCP_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==3));
    pRP_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==3));
    pSC_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==1 & FAgmstd_grp==3));
    pSR_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==3 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==3));
    pSP_RESOLVE = ranksum(FAgmstd_boxplot(subject_grp==4 & FAgmstd_grp==3),FAgmstd_boxplot(subject_grp==2 & FAgmstd_grp==3));
    table_pvals(10,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(11).fig = figure(11);
    set(h(11).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAgmstd_grp)/2-0.20 unique(FAgmstd_grp)/2+0.20]', repmat(FAgmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAgmstd_grp)/2-0.20 unique(FAgmstd_grp)/2+0.20]', repmat(FAgmstd_median',1,2)','m-','LineWidth',8)
    scatter(FAgmstd_grp(subject_grp==3)/2, FAgmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAgmstd_grp(subject_grp==1)/2, FAgmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAgmstd_grp(subject_grp==4)/2, FAgmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAgmstd_grp(subject_grp==2)/2, FAgmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.185,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.180,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.175,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.170,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.165,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.160,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.185,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.180,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.175,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.170,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.165,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.160,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.185,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.185,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.180,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.180,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.175,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.170,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.170,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.165,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.165,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.160,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA standard deviation value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.06 0.190])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_std_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEAN VALUE DISTRIBUTIONS FROM C3-C6
    pCR_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==1),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==2),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==1 & fieldcoefmean_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==3 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmean_boxplot(subject_grp==4 & fieldcoefmean_grp==3),fieldcoefmean_boxplot(subject_grp==2 & fieldcoefmean_grp==3));
    table_pvals(11,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(12).fig = figure(12);
    set(h(12).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmean_grp)/2-0.20 unique(fieldcoefmean_grp)/2+0.20]', repmat(fieldcoefmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmean_grp)/2-0.20 unique(fieldcoefmean_grp)/2+0.20]', repmat(fieldcoefmean_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmean_grp(subject_grp==3)/2, fieldcoefmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean_grp(subject_grp==1)/2, fieldcoefmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean_grp(subject_grp==4)/2, fieldcoefmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean_grp(subject_grp==2)/2, fieldcoefmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| mean value distribution from';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| mean values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_mean_C3C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEDIAN VALUE DISTRIBUTIONS FROM C3-C6
    pCR_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==1),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==2),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==1 & fieldcoefmedian_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==3 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmedian_boxplot(subject_grp==4 & fieldcoefmedian_grp==3),fieldcoefmedian_boxplot(subject_grp==2 & fieldcoefmedian_grp==3));
    table_pvals(12,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(13).fig = figure(13);
    set(h(13).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmedian_grp)/2-0.20 unique(fieldcoefmedian_grp)/2+0.20]', repmat(fieldcoefmedian_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmedian_grp)/2-0.20 unique(fieldcoefmedian_grp)/2+0.20]', repmat(fieldcoefmedian_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmedian_grp(subject_grp==3)/2, fieldcoefmedian_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian_grp(subject_grp==1)/2, fieldcoefmedian_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian_grp(subject_grp==4)/2, fieldcoefmedian_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian_grp(subject_grp==2)/2, fieldcoefmedian_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| median value distribution from';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| median values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_median_C3C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF STD VALUE DISTRIBUTIONS FROM C3-C6
    pCR_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==1),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==2),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==1 & fieldcoefstd_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==3 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefstd_boxplot(subject_grp==4 & fieldcoefstd_grp==3),fieldcoefstd_boxplot(subject_grp==2 & fieldcoefstd_grp==3));
    table_pvals(13,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(14).fig = figure(14);
    set(h(14).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefstd_grp)/2-0.20 unique(fieldcoefstd_grp)/2+0.20]', repmat(fieldcoefstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefstd_grp)/2-0.20 unique(fieldcoefstd_grp)/2+0.20]', repmat(fieldcoefstd_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefstd_grp(subject_grp==3)/2, fieldcoefstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd_grp(subject_grp==1)/2, fieldcoefstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd_grp(subject_grp==4)/2, fieldcoefstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd_grp(subject_grp==2)/2, fieldcoefstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| standard deviation value';'distribution from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| SD values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_std_C3C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEAN VALUE DISTRIBUTIONS FROM C3
    pCR_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==1),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==2),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==1 & fieldcoefmean3_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==3 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmean3_boxplot(subject_grp==4 & fieldcoefmean3_grp==3),fieldcoefmean3_boxplot(subject_grp==2 & fieldcoefmean3_grp==3));
    table_pvals(14,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(15).fig = figure(15);
    set(h(15).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmean3_grp)/2-0.20 unique(fieldcoefmean3_grp)/2+0.20]', repmat(fieldcoefmean3_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmean3_grp)/2-0.20 unique(fieldcoefmean3_grp)/2+0.20]', repmat(fieldcoefmean3_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmean3_grp(subject_grp==3)/2, fieldcoefmean3_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean3_grp(subject_grp==1)/2, fieldcoefmean3_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean3_grp(subject_grp==4)/2, fieldcoefmean3_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean3_grp(subject_grp==2)/2, fieldcoefmean3_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| mean value distribution from';'C3 segment over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| mean values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_mean_C3_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEDIAN VALUE DISTRIBUTIONS FROM C3
    pCR_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==1),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==2),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==1 & fieldcoefmedian3_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==3 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmedian3_boxplot(subject_grp==4 & fieldcoefmedian3_grp==3),fieldcoefmedian3_boxplot(subject_grp==2 & fieldcoefmedian3_grp==3));
    table_pvals(15,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(16).fig = figure(16);
    set(h(16).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmedian3_grp)/2-0.20 unique(fieldcoefmedian3_grp)/2+0.20]', repmat(fieldcoefmedian3_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmedian3_grp)/2-0.20 unique(fieldcoefmedian3_grp)/2+0.20]', repmat(fieldcoefmedian3_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmedian3_grp(subject_grp==3)/2, fieldcoefmedian3_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian3_grp(subject_grp==1)/2, fieldcoefmedian3_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian3_grp(subject_grp==4)/2, fieldcoefmedian3_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian3_grp(subject_grp==2)/2, fieldcoefmedian3_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| median value distribution from';'C3 segment over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| median values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_median_C3_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF STD VALUE DISTRIBUTIONS FROM C3
    pCR_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==1),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==2),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==1 & fieldcoefstd3_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==3 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefstd3_boxplot(subject_grp==4 & fieldcoefstd3_grp==3),fieldcoefstd3_boxplot(subject_grp==2 & fieldcoefstd3_grp==3));
    table_pvals(16,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(17).fig = figure(17);
    set(h(17).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefstd3_grp)/2-0.20 unique(fieldcoefstd3_grp)/2+0.20]', repmat(fieldcoefstd3_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefstd3_grp)/2-0.20 unique(fieldcoefstd3_grp)/2+0.20]', repmat(fieldcoefstd3_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefstd3_grp(subject_grp==3)/2, fieldcoefstd3_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd3_grp(subject_grp==1)/2, fieldcoefstd3_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd3_grp(subject_grp==4)/2, fieldcoefstd3_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd3_grp(subject_grp==2)/2, fieldcoefstd3_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| standard deviation value';'distribution from C3 segment over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| standard deviation values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_std_C3_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEAN VALUE DISTRIBUTIONS FROM C5-C6
    pCR_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==1),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==2),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==1 & fieldcoefmean56_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==3 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmean56_boxplot(subject_grp==4 & fieldcoefmean56_grp==3),fieldcoefmean56_boxplot(subject_grp==2 & fieldcoefmean56_grp==3));
    table_pvals(17,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(18).fig = figure(18);
    set(h(18).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmean56_grp)/2-0.20 unique(fieldcoefmean56_grp)/2+0.20]', repmat(fieldcoefmean56_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmean56_grp)/2-0.20 unique(fieldcoefmean56_grp)/2+0.20]', repmat(fieldcoefmean56_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmean56_grp(subject_grp==3)/2, fieldcoefmean56_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean56_grp(subject_grp==1)/2, fieldcoefmean56_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean56_grp(subject_grp==4)/2, fieldcoefmean56_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmean56_grp(subject_grp==2)/2, fieldcoefmean56_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| mean value distribution from';'C5-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| mean values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_mean_C5C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF MEDIAN VALUE DISTRIBUTIONS FROM C5-C6
    pCR_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==1),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==2),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==1 & fieldcoefmedian56_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==3 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefmedian56_boxplot(subject_grp==4 & fieldcoefmedian56_grp==3),fieldcoefmedian56_boxplot(subject_grp==2 & fieldcoefmedian56_grp==3));
    table_pvals(18,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(19).fig = figure(19);
    set(h(19).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefmedian56_grp)/2-0.20 unique(fieldcoefmedian56_grp)/2+0.20]', repmat(fieldcoefmedian56_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefmedian56_grp)/2-0.20 unique(fieldcoefmedian56_grp)/2+0.20]', repmat(fieldcoefmedian56_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefmedian56_grp(subject_grp==3)/2, fieldcoefmedian56_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian56_grp(subject_grp==1)/2, fieldcoefmedian56_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian56_grp(subject_grp==4)/2, fieldcoefmedian56_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefmedian56_grp(subject_grp==2)/2, fieldcoefmedian56_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| median value distribution from';'C5-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| median values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_median_C5C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FIELDCOEF STD VALUE DISTRIBUTIONS FROM C5-C6
    pCR_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==1));
    pCP_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==1));
    pRP_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==1));
    pSC_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==1));
    pSR_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==1));
    pSP_ZOOMint = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==1),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==1));
    pCR_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==2));
    pCP_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==2));
    pRP_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==2));
    pSC_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==2));
    pSR_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==2));
    pSP_ZOOMnotint = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==2),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==2));
    pCR_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==3));
    pCP_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==3));
    pRP_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==3));
    pSC_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==1 & fieldcoefstd56_grp==3));
    pSR_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==3 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==3));
    pSP_RESOLVE = ranksum(fieldcoefstd56_boxplot(subject_grp==4 & fieldcoefstd56_grp==3),fieldcoefstd56_boxplot(subject_grp==2 & fieldcoefstd56_grp==3));
    table_pvals(19,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(20).fig = figure(20);
    set(h(20).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(fieldcoefstd56_grp)/2-0.20 unique(fieldcoefstd56_grp)/2+0.20]', repmat(fieldcoefstd56_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(fieldcoefstd56_grp)/2-0.20 unique(fieldcoefstd56_grp)/2+0.20]', repmat(fieldcoefstd56_median',1,2)','m-','LineWidth',8)
    scatter(fieldcoefstd56_grp(subject_grp==3)/2, fieldcoefstd56_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd56_grp(subject_grp==1)/2, fieldcoefstd56_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd56_grp(subject_grp==4)/2, fieldcoefstd56_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(fieldcoefstd56_grp(subject_grp==2)/2, fieldcoefstd56_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,79.0,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,77.0,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,75.0,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,73.0,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,71.0,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,69.0,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,79.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,79.0,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,77.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,77.0,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,75.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,75.0,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,73.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,73.0,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,71.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,71.0,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,69.0,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,69.0,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Topup-output |fieldcoef| standard deviation value';'distribution from C5-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject |fieldcoef| standard deviation values [Hz]','FontSize',18)
    axis([0.5/2 3.5/2 0 80])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_fieldcoef_std_C5C6_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEDIAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==1),MDWM_boxplot(subject_grp==3 & MDWM_grp==1));
    pCP_ZOOMint = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==1),MDWM_boxplot(subject_grp==2 & MDWM_grp==1));
    pRP_ZOOMint = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==1),MDWM_boxplot(subject_grp==2 & MDWM_grp==1));
    pSC_ZOOMint = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==1),MDWM_boxplot(subject_grp==1 & MDWM_grp==1));
    pSR_ZOOMint = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==1),MDWM_boxplot(subject_grp==4 & MDWM_grp==1));
    pSP_ZOOMint = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==1),MDWM_boxplot(subject_grp==2 & MDWM_grp==1));
    pCR_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==2),MDWM_boxplot(subject_grp==3 & MDWM_grp==2));
    pCP_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==2),MDWM_boxplot(subject_grp==2 & MDWM_grp==2));
    pRP_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==2),MDWM_boxplot(subject_grp==2 & MDWM_grp==2));
    pSC_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==2),MDWM_boxplot(subject_grp==1 & MDWM_grp==2));
    pSR_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==2),MDWM_boxplot(subject_grp==4 & MDWM_grp==2));
    pSP_ZOOMnotint = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==2),MDWM_boxplot(subject_grp==2 & MDWM_grp==2));
    pCR_RESOLVE = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==3),MDWM_boxplot(subject_grp==3 & MDWM_grp==3));
    pCP_RESOLVE = ranksum(MDWM_boxplot(subject_grp==1 & MDWM_grp==3),MDWM_boxplot(subject_grp==2 & MDWM_grp==3));
    pRP_RESOLVE = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==3),MDWM_boxplot(subject_grp==2 & MDWM_grp==3));
    pSC_RESOLVE = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==3),MDWM_boxplot(subject_grp==1 & MDWM_grp==3));
    pSR_RESOLVE = ranksum(MDWM_boxplot(subject_grp==3 & MDWM_grp==3),MDWM_boxplot(subject_grp==4 & MDWM_grp==3));
    pSP_RESOLVE = ranksum(MDWM_boxplot(subject_grp==4 & MDWM_grp==3),MDWM_boxplot(subject_grp==2 & MDWM_grp==3));
    table_pvals(20,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(25).fig = figure(25);
    set(h(25).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWM_grp)/2-0.20 unique(MDWM_grp)/2+0.20]', repmat(MDWM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWM_grp)/2-0.20 unique(MDWM_grp)/2+0.20]', repmat(MDWM_median',1,2)','m-','LineWidth',8)
    scatter(MDWM_grp(subject_grp==3)/2, MDWM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWM_grp(subject_grp==1)/2, MDWM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWM_grp(subject_grp==4)/2, MDWM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWM_grp(subject_grp==2)/2, MDWM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD median value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD median values','FontSize',18)
    axis([0.5/2 3.5/2 0.7 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_median_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEDIAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==1),MDGM_boxplot(subject_grp==3 & MDGM_grp==1));
    pCP_ZOOMint = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==1),MDGM_boxplot(subject_grp==2 & MDGM_grp==1));
    pRP_ZOOMint = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==1),MDGM_boxplot(subject_grp==2 & MDGM_grp==1));
    pSC_ZOOMint = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==1),MDGM_boxplot(subject_grp==1 & MDGM_grp==1));
    pSR_ZOOMint = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==1),MDGM_boxplot(subject_grp==4 & MDGM_grp==1));
    pSP_ZOOMint = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==1),MDGM_boxplot(subject_grp==2 & MDGM_grp==1));
    pCR_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==2),MDGM_boxplot(subject_grp==3 & MDGM_grp==2));
    pCP_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==2),MDGM_boxplot(subject_grp==2 & MDGM_grp==2));
    pRP_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==2),MDGM_boxplot(subject_grp==2 & MDGM_grp==2));
    pSC_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==2),MDGM_boxplot(subject_grp==1 & MDGM_grp==2));
    pSR_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==2),MDGM_boxplot(subject_grp==4 & MDGM_grp==2));
    pSP_ZOOMnotint = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==2),MDGM_boxplot(subject_grp==2 & MDGM_grp==2));
    pCR_RESOLVE = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==3),MDGM_boxplot(subject_grp==3 & MDGM_grp==3));
    pCP_RESOLVE = ranksum(MDGM_boxplot(subject_grp==1 & MDGM_grp==3),MDGM_boxplot(subject_grp==2 & MDGM_grp==3));
    pRP_RESOLVE = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==3),MDGM_boxplot(subject_grp==2 & MDGM_grp==3));
    pSC_RESOLVE = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==3),MDGM_boxplot(subject_grp==1 & MDGM_grp==3));
    pSR_RESOLVE = ranksum(MDGM_boxplot(subject_grp==3 & MDGM_grp==3),MDGM_boxplot(subject_grp==4 & MDGM_grp==3));
    pSP_RESOLVE = ranksum(MDGM_boxplot(subject_grp==4 & MDGM_grp==3),MDGM_boxplot(subject_grp==2 & MDGM_grp==3));
    table_pvals(21,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(26).fig = figure(26);
    set(h(26).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDGM_grp)/2-0.20 unique(MDGM_grp)/2+0.20]', repmat(MDGM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDGM_grp)/2-0.20 unique(MDGM_grp)/2+0.20]', repmat(MDGM_median',1,2)','m-','LineWidth',8)
    scatter(MDGM_grp(subject_grp==3)/2, MDGM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGM_grp(subject_grp==1)/2, MDGM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGM_grp(subject_grp==4)/2, MDGM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDGM_grp(subject_grp==2)/2, MDGM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.550,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.450,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD median value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD median values','FontSize',18)
    axis([0.5/2 3.5/2 0.7 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_median_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEDIAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==1));
    pCP_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==1));
    pRP_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==1));
    pSC_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==1));
    pSR_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==1));
    pSP_ZOOMint = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==1),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==2),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==2));
    pCR_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==3));
    pCP_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==3));
    pRP_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==3));
    pSC_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==1 & MDWMGMdiff_grp==3));
    pSR_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==3 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==3));
    pSP_RESOLVE = ranksum(MDWMGMdiff_boxplot(subject_grp==4 & MDWMGMdiff_grp==3),MDWMGMdiff_boxplot(subject_grp==2 & MDWMGMdiff_grp==3));
    table_pvals(22,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(27).fig = figure(27);
    set(h(27).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMGMdiff_grp)/2-0.20 unique(MDWMGMdiff_grp)/2+0.20]', repmat(MDWMGMdiff_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMGMdiff_grp)/2-0.20 unique(MDWMGMdiff_grp)/2+0.20]', repmat(MDWMGMdiff_median',1,2)','m-','LineWidth',8)
    scatter(MDWMGMdiff_grp(subject_grp==3)/2, MDWMGMdiff_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiff_grp(subject_grp==1)/2, MDWMGMdiff_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiff_grp(subject_grp==4)/2, MDWMGMdiff_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiff_grp(subject_grp==2)/2, MDWMGMdiff_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMint*1000000)/1000000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.49,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.47,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.45,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.43,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.41,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.39,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'median MD values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD median diferences','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.50])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_median_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD STD VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==1));
    pCP_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==1));
    pRP_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==1));
    pSC_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==1));
    pSR_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==1));
    pSP_ZOOMint = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==1),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==1));
    pCR_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==2));
    pCP_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==2));
    pRP_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==2));
    pSC_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==2));
    pSR_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==2));
    pSP_ZOOMnotint = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==2),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==2));
    pCR_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==3));
    pCP_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==3));
    pRP_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==3));
    pSC_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==1 & MDwmstd_grp==3));
    pSR_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==3 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==3));
    pSP_RESOLVE = ranksum(MDwmstd_boxplot(subject_grp==4 & MDwmstd_grp==3),MDwmstd_boxplot(subject_grp==2 & MDwmstd_grp==3));
    table_pvals(23,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(28).fig = figure(28);
    set(h(28).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDwmstd_grp)/2-0.20 unique(MDwmstd_grp)/2+0.20]', repmat(MDwmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDwmstd_grp)/2-0.20 unique(MDwmstd_grp)/2+0.20]', repmat(MDwmstd_median',1,2)','m-','LineWidth',8)
    scatter(MDwmstd_grp(subject_grp==3)/2, MDwmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDwmstd_grp(subject_grp==1)/2, MDwmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDwmstd_grp(subject_grp==4)/2, MDwmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDwmstd_grp(subject_grp==2)/2, MDwmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.790,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.760,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.670,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.640,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD standard deviation value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.00 0.800])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_std_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD STD VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==1));
    pCP_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==1));
    pRP_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==1));
    pSC_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==1));
    pSR_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==1));
    pSP_ZOOMint = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==1),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==1));
    pCR_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==2));
    pCP_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==2));
    pRP_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==2));
    pSC_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==2));
    pSR_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==2));
    pSP_ZOOMnotint = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==2),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==2));
    pCR_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==3));
    pCP_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==3));
    pRP_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==3));
    pSC_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==1 & MDgmstd_grp==3));
    pSR_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==3 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==3));
    pSP_RESOLVE = ranksum(MDgmstd_boxplot(subject_grp==4 & MDgmstd_grp==3),MDgmstd_boxplot(subject_grp==2 & MDgmstd_grp==3));
    table_pvals(24,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(29).fig = figure(29);
    set(h(29).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDgmstd_grp)/2-0.20 unique(MDgmstd_grp)/2+0.20]', repmat(MDgmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDgmstd_grp)/2-0.20 unique(MDgmstd_grp)/2+0.20]', repmat(MDgmstd_median',1,2)','m-','LineWidth',8)
    scatter(MDgmstd_grp(subject_grp==3)/2, MDgmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDgmstd_grp(subject_grp==1)/2, MDgmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDgmstd_grp(subject_grp==4)/2, MDgmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDgmstd_grp(subject_grp==2)/2, MDgmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.790,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.760,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.670,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.640,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD standard deviation value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.00 0.800])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_std_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEDIAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==1),f1WM_boxplot(subject_grp==3 & f1WM_grp==1));
    pCP_ZOOMint = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==1),f1WM_boxplot(subject_grp==2 & f1WM_grp==1));
    pRP_ZOOMint = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==1),f1WM_boxplot(subject_grp==2 & f1WM_grp==1));
    pSC_ZOOMint = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==1),f1WM_boxplot(subject_grp==1 & f1WM_grp==1));
    pSR_ZOOMint = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==1),f1WM_boxplot(subject_grp==4 & f1WM_grp==1));
    pSP_ZOOMint = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==1),f1WM_boxplot(subject_grp==2 & f1WM_grp==1));
    pCR_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==2),f1WM_boxplot(subject_grp==3 & f1WM_grp==2));
    pCP_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==2),f1WM_boxplot(subject_grp==2 & f1WM_grp==2));
    pRP_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==2),f1WM_boxplot(subject_grp==2 & f1WM_grp==2));
    pSC_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==2),f1WM_boxplot(subject_grp==1 & f1WM_grp==2));
    pSR_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==2),f1WM_boxplot(subject_grp==4 & f1WM_grp==2));
    pSP_ZOOMnotint = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==2),f1WM_boxplot(subject_grp==2 & f1WM_grp==2));
    pCR_RESOLVE = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==3),f1WM_boxplot(subject_grp==3 & f1WM_grp==3));
    pCP_RESOLVE = ranksum(f1WM_boxplot(subject_grp==1 & f1WM_grp==3),f1WM_boxplot(subject_grp==2 & f1WM_grp==3));
    pRP_RESOLVE = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==3),f1WM_boxplot(subject_grp==2 & f1WM_grp==3));
    pSC_RESOLVE = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==3),f1WM_boxplot(subject_grp==1 & f1WM_grp==3));
    pSR_RESOLVE = ranksum(f1WM_boxplot(subject_grp==3 & f1WM_grp==3),f1WM_boxplot(subject_grp==4 & f1WM_grp==3));
    pSP_RESOLVE = ranksum(f1WM_boxplot(subject_grp==4 & f1WM_grp==3),f1WM_boxplot(subject_grp==2 & f1WM_grp==3));
    table_pvals(25,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(30).fig = figure(30);
    set(h(30).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WM_grp)/2-0.20 unique(f1WM_grp)/2+0.20]', repmat(f1WM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WM_grp)/2-0.20 unique(f1WM_grp)/2+0.20]', repmat(f1WM_median',1,2)','m-','LineWidth',8)
    scatter(f1WM_grp(subject_grp==3)/2, f1WM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WM_grp(subject_grp==1)/2, f1WM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WM_grp(subject_grp==4)/2, f1WM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WM_grp(subject_grp==2)/2, f1WM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.680,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.655,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.630,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.605,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.580,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.555,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 median value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 median values','FontSize',18)
    axis([0.5/2 3.5/2 0.25 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_median_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEDIAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==1),f1GM_boxplot(subject_grp==3 & f1GM_grp==1));
    pCP_ZOOMint = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==1),f1GM_boxplot(subject_grp==2 & f1GM_grp==1));
    pRP_ZOOMint = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==1),f1GM_boxplot(subject_grp==2 & f1GM_grp==1));
    pSC_ZOOMint = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==1),f1GM_boxplot(subject_grp==1 & f1GM_grp==1));
    pSR_ZOOMint = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==1),f1GM_boxplot(subject_grp==4 & f1GM_grp==1));
    pSP_ZOOMint = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==1),f1GM_boxplot(subject_grp==2 & f1GM_grp==1));
    pCR_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==2),f1GM_boxplot(subject_grp==3 & f1GM_grp==2));
    pCP_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==2),f1GM_boxplot(subject_grp==2 & f1GM_grp==2));
    pRP_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==2),f1GM_boxplot(subject_grp==2 & f1GM_grp==2));
    pSC_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==2),f1GM_boxplot(subject_grp==1 & f1GM_grp==2));
    pSR_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==2),f1GM_boxplot(subject_grp==4 & f1GM_grp==2));
    pSP_ZOOMnotint = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==2),f1GM_boxplot(subject_grp==2 & f1GM_grp==2));
    pCR_RESOLVE = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==3),f1GM_boxplot(subject_grp==3 & f1GM_grp==3));
    pCP_RESOLVE = ranksum(f1GM_boxplot(subject_grp==1 & f1GM_grp==3),f1GM_boxplot(subject_grp==2 & f1GM_grp==3));
    pRP_RESOLVE = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==3),f1GM_boxplot(subject_grp==2 & f1GM_grp==3));
    pSC_RESOLVE = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==3),f1GM_boxplot(subject_grp==1 & f1GM_grp==3));
    pSR_RESOLVE = ranksum(f1GM_boxplot(subject_grp==3 & f1GM_grp==3),f1GM_boxplot(subject_grp==4 & f1GM_grp==3));
    pSP_RESOLVE = ranksum(f1GM_boxplot(subject_grp==4 & f1GM_grp==3),f1GM_boxplot(subject_grp==2 & f1GM_grp==3));
    table_pvals(26,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(31).fig = figure(31);
    set(h(31).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1GM_grp)/2-0.20 unique(f1GM_grp)/2+0.20]', repmat(f1GM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1GM_grp)/2-0.20 unique(f1GM_grp)/2+0.20]', repmat(f1GM_median',1,2)','m-','LineWidth',8)
    scatter(f1GM_grp(subject_grp==3)/2, f1GM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GM_grp(subject_grp==1)/2, f1GM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GM_grp(subject_grp==4)/2, f1GM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1GM_grp(subject_grp==2)/2, f1GM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.680,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.655,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.630,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.605,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.580,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.555,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 median value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 median values','FontSize',18)
    axis([0.5/2 3.5/2 0.25 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_median_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEDIAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==1));
    pCP_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==1));
    pRP_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==1));
    pSC_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==1));
    pSR_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==1));
    pSP_ZOOMint = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==1),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==2),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==2));
    pCR_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==3));
    pCP_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==3));
    pRP_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==3));
    pSC_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==1 & f1WMGMdiff_grp==3));
    pSR_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==3 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==3));
    pSP_RESOLVE = ranksum(f1WMGMdiff_boxplot(subject_grp==4 & f1WMGMdiff_grp==3),f1WMGMdiff_boxplot(subject_grp==2 & f1WMGMdiff_grp==3));
    table_pvals(27,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(32).fig = figure(32);
    set(h(32).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMGMdiff_grp)/2-0.20 unique(f1WMGMdiff_grp)/2+0.20]', repmat(f1WMGMdiff_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMGMdiff_grp)/2-0.20 unique(f1WMGMdiff_grp)/2+0.20]', repmat(f1WMGMdiff_median',1,2)','m-','LineWidth',8)
    scatter(f1WMGMdiff_grp(subject_grp==3)/2, f1WMGMdiff_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiff_grp(subject_grp==1)/2, f1WMGMdiff_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiff_grp(subject_grp==4)/2, f1WMGMdiff_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiff_grp(subject_grp==2)/2, f1WMGMdiff_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.190,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.190,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.175,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.160,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.145,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.145,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.130,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.130,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.115,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.115,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.190,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.190,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.175,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.160,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.145,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.145,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.130,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.130,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.115,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.115,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.190,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.190,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.175,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.175,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.160,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.160,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.145,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.145,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.130,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.130,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.115,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.115,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'median f_1 values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject f_1 median values','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_median_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 STD VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==1));
    pCP_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==1));
    pRP_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==1));
    pSC_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==1));
    pSR_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==1));
    pSP_ZOOMint = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==1),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==1));
    pCR_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==2));
    pCP_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==2));
    pRP_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==2));
    pSC_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==2));
    pSR_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==2));
    pSP_ZOOMnotint = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==2),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==2));
    pCR_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==3));
    pCP_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==3));
    pRP_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==3));
    pSC_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==1 & f1wmstd_grp==3));
    pSR_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==3 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==3));
    pSP_RESOLVE = ranksum(f1wmstd_boxplot(subject_grp==4 & f1wmstd_grp==3),f1wmstd_boxplot(subject_grp==2 & f1wmstd_grp==3));
    table_pvals(28,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(33).fig = figure(33);
    set(h(33).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1wmstd_grp)/2-0.20 unique(f1wmstd_grp)/2+0.20]', repmat(f1wmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1wmstd_grp)/2-0.20 unique(f1wmstd_grp)/2+0.20]', repmat(f1wmstd_median',1,2)','m-','LineWidth',8)
    scatter(f1wmstd_grp(subject_grp==3)/2, f1wmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1wmstd_grp(subject_grp==1)/2, f1wmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1wmstd_grp(subject_grp==4)/2, f1wmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1wmstd_grp(subject_grp==2)/2, f1wmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.217,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.212,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.205,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.205,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.198,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.191,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.191,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.193,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.217,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.212,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.207,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.207,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.202,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.202,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.198,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.193,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.217,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.212,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.207,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.207,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.202,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.202,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.198,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.193,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 standard deviation value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.06 0.220])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_std_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 STD VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==1));
    pCP_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==1));
    pRP_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==1));
    pSC_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==1));
    pSR_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==1));
    pSP_ZOOMint = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==1),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==1));
    pCR_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==2));
    pCP_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==2));
    pRP_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==2));
    pSC_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==2));
    pSR_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==2));
    pSP_ZOOMnotint = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==2),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==2));
    pCR_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==3));
    pCP_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==3));
    pRP_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==3));
    pSC_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==1 & f1gmstd_grp==3));
    pSR_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==3 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==3));
    pSP_RESOLVE = ranksum(f1gmstd_boxplot(subject_grp==4 & f1gmstd_grp==3),f1gmstd_boxplot(subject_grp==2 & f1gmstd_grp==3));
    table_pvals(29,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(34).fig = figure(34);
    set(h(34).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1gmstd_grp)/2-0.20 unique(f1gmstd_grp)/2+0.20]', repmat(f1gmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1gmstd_grp)/2-0.20 unique(f1gmstd_grp)/2+0.20]', repmat(f1gmstd_median',1,2)','m-','LineWidth',8)
    scatter(f1gmstd_grp(subject_grp==3)/2, f1gmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1gmstd_grp(subject_grp==1)/2, f1gmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1gmstd_grp(subject_grp==4)/2, f1gmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1gmstd_grp(subject_grp==2)/2, f1gmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.217,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.212,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.207,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.207,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.202,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.202,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.198,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.193,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.217,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.212,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.207,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.207,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.202,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.202,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.198,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.193,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.217,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.217,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.212,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.212,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.207,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.207,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.202,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.202,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.198,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.198,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.193,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.193,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 standard deviation value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.06 0.220])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_std_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEDIAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==1),dWM_boxplot(subject_grp==3 & dWM_grp==1));
    pCP_ZOOMint = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==1),dWM_boxplot(subject_grp==2 & dWM_grp==1));
    pRP_ZOOMint = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==1),dWM_boxplot(subject_grp==2 & dWM_grp==1));
    pSC_ZOOMint = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==1),dWM_boxplot(subject_grp==1 & dWM_grp==1));
    pSR_ZOOMint = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==1),dWM_boxplot(subject_grp==4 & dWM_grp==1));
    pSP_ZOOMint = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==1),dWM_boxplot(subject_grp==2 & dWM_grp==1));
    pCR_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==2),dWM_boxplot(subject_grp==3 & dWM_grp==2));
    pCP_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==2),dWM_boxplot(subject_grp==2 & dWM_grp==2));
    pRP_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==2),dWM_boxplot(subject_grp==2 & dWM_grp==2));
    pSC_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==2),dWM_boxplot(subject_grp==1 & dWM_grp==2));
    pSR_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==2),dWM_boxplot(subject_grp==4 & dWM_grp==2));
    pSP_ZOOMnotint = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==2),dWM_boxplot(subject_grp==2 & dWM_grp==2));
    pCR_RESOLVE = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==3),dWM_boxplot(subject_grp==3 & dWM_grp==3));
    pCP_RESOLVE = ranksum(dWM_boxplot(subject_grp==1 & dWM_grp==3),dWM_boxplot(subject_grp==2 & dWM_grp==3));
    pRP_RESOLVE = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==3),dWM_boxplot(subject_grp==2 & dWM_grp==3));
    pSC_RESOLVE = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==3),dWM_boxplot(subject_grp==1 & dWM_grp==3));
    pSR_RESOLVE = ranksum(dWM_boxplot(subject_grp==3 & dWM_grp==3),dWM_boxplot(subject_grp==4 & dWM_grp==3));
    pSP_RESOLVE = ranksum(dWM_boxplot(subject_grp==4 & dWM_grp==3),dWM_boxplot(subject_grp==2 & dWM_grp==3));
    table_pvals(30,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(35).fig = figure(35);
    set(h(35).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWM_grp)/2-0.20 unique(dWM_grp)/2+0.20]', repmat(dWM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWM_grp)/2-0.20 unique(dWM_grp)/2+0.20]', repmat(dWM_median',1,2)','m-','LineWidth',8)
    scatter(dWM_grp(subject_grp==3)/2, dWM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWM_grp(subject_grp==1)/2, dWM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWM_grp(subject_grp==4)/2, dWM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWM_grp(subject_grp==2)/2, dWM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d median value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d median values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_median_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEDIAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==1),dGM_boxplot(subject_grp==3 & dGM_grp==1));
    pCP_ZOOMint = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==1),dGM_boxplot(subject_grp==2 & dGM_grp==1));
    pRP_ZOOMint = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==1),dGM_boxplot(subject_grp==2 & dGM_grp==1));
    pSC_ZOOMint = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==1),dGM_boxplot(subject_grp==1 & dGM_grp==1));
    pSR_ZOOMint = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==1),dGM_boxplot(subject_grp==4 & dGM_grp==1));
    pSP_ZOOMint = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==1),dGM_boxplot(subject_grp==2 & dGM_grp==1));
    pCR_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==2),dGM_boxplot(subject_grp==3 & dGM_grp==2));
    pCP_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==2),dGM_boxplot(subject_grp==2 & dGM_grp==2));
    pRP_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==2),dGM_boxplot(subject_grp==2 & dGM_grp==2));
    pSC_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==2),dGM_boxplot(subject_grp==1 & dGM_grp==2));
    pSR_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==2),dGM_boxplot(subject_grp==4 & dGM_grp==2));
    pSP_ZOOMnotint = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==2),dGM_boxplot(subject_grp==2 & dGM_grp==2));
    pCR_RESOLVE = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==3),dGM_boxplot(subject_grp==3 & dGM_grp==3));
    pCP_RESOLVE = ranksum(dGM_boxplot(subject_grp==1 & dGM_grp==3),dGM_boxplot(subject_grp==2 & dGM_grp==3));
    pRP_RESOLVE = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==3),dGM_boxplot(subject_grp==2 & dGM_grp==3));
    pSC_RESOLVE = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==3),dGM_boxplot(subject_grp==1 & dGM_grp==3));
    pSR_RESOLVE = ranksum(dGM_boxplot(subject_grp==3 & dGM_grp==3),dGM_boxplot(subject_grp==4 & dGM_grp==3));
    pSP_RESOLVE = ranksum(dGM_boxplot(subject_grp==4 & dGM_grp==3),dGM_boxplot(subject_grp==2 & dGM_grp==3));
    table_pvals(31,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(36).fig = figure(36);
    set(h(36).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dGM_grp)/2-0.20 unique(dGM_grp)/2+0.20]', repmat(dGM_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dGM_grp)/2-0.20 unique(dGM_grp)/2+0.20]', repmat(dGM_median',1,2)','m-','LineWidth',8)
    scatter(dGM_grp(subject_grp==3)/2, dGM_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGM_grp(subject_grp==1)/2, dGM_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGM_grp(subject_grp==4)/2, dGM_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dGM_grp(subject_grp==2)/2, dGM_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d median value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d median values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_median_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEDIAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==1));
    pCP_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==1));
    pRP_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==1));
    pSC_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==1));
    pSR_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==1));
    pSP_ZOOMint = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==1),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==1));
    pCR_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==2));
    pCP_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==2));
    pRP_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==2));
    pSC_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==2));
    pSR_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==2));
    pSP_ZOOMnotint = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==2),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==2));
    pCR_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==3));
    pCP_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==3));
    pRP_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==3));
    pSC_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==1 & dWMGMdiff_grp==3));
    pSR_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==3 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==3));
    pSP_RESOLVE = ranksum(dWMGMdiff_boxplot(subject_grp==4 & dWMGMdiff_grp==3),dWMGMdiff_boxplot(subject_grp==2 & dWMGMdiff_grp==3));
    table_pvals(32,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(37).fig = figure(37);
    set(h(37).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMGMdiff_grp)/2-0.20 unique(dWMGMdiff_grp)/2+0.20]', repmat(dWMGMdiff_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMGMdiff_grp)/2-0.20 unique(dWMGMdiff_grp)/2+0.20]', repmat(dWMGMdiff_median',1,2)','m-','LineWidth',8)
    scatter(dWMGMdiff_grp(subject_grp==3)/2, dWMGMdiff_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiff_grp(subject_grp==1)/2, dWMGMdiff_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiff_grp(subject_grp==4)/2, dWMGMdiff_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiff_grp(subject_grp==2)/2, dWMGMdiff_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.83,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.80,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.77,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.74,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.71,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.68,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'median d values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d median differences','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.85])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_median_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d STD VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==3 & dwmstd_grp==1));
    pCP_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==1));
    pRP_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==1));
    pSC_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==1 & dwmstd_grp==1));
    pSR_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==4 & dwmstd_grp==1));
    pSP_ZOOMint = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==1),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==1));
    pCR_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==3 & dwmstd_grp==2));
    pCP_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==2));
    pRP_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==2));
    pSC_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==1 & dwmstd_grp==2));
    pSR_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==4 & dwmstd_grp==2));
    pSP_ZOOMnotint = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==2),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==2));
    pCR_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==3 & dwmstd_grp==3));
    pCP_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==1 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==3));
    pRP_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==3));
    pSC_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==1 & dwmstd_grp==3));
    pSR_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==3 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==4 & dwmstd_grp==3));
    pSP_RESOLVE = ranksum(dwmstd_boxplot(subject_grp==4 & dwmstd_grp==3),dwmstd_boxplot(subject_grp==2 & dwmstd_grp==3));
    table_pvals(33,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(38).fig = figure(38);
    set(h(38).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dwmstd_grp)/2-0.20 unique(dwmstd_grp)/2+0.20]', repmat(dwmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dwmstd_grp)/2-0.20 unique(dwmstd_grp)/2+0.20]', repmat(dwmstd_median',1,2)','m-','LineWidth',8)
    scatter(dwmstd_grp(subject_grp==3)/2, dwmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dwmstd_grp(subject_grp==1)/2, dwmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dwmstd_grp(subject_grp==4)/2, dwmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dwmstd_grp(subject_grp==2)/2, dwmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.790,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.760,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.670,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.640,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d standard deviation value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.00 0.800])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_std_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d STD VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==3 & dgmstd_grp==1));
    pCP_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==1));
    pRP_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==1));
    pSC_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==1 & dgmstd_grp==1));
    pSR_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==4 & dgmstd_grp==1));
    pSP_ZOOMint = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==1),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==1));
    pCR_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==3 & dgmstd_grp==2));
    pCP_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==2));
    pRP_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==2));
    pSC_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==1 & dgmstd_grp==2));
    pSR_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==4 & dgmstd_grp==2));
    pSP_ZOOMnotint = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==2),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==2));
    pCR_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==3 & dgmstd_grp==3));
    pCP_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==1 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==3));
    pRP_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==3));
    pSC_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==1 & dgmstd_grp==3));
    pSR_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==3 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==4 & dgmstd_grp==3));
    pSP_RESOLVE = ranksum(dgmstd_boxplot(subject_grp==4 & dgmstd_grp==3),dgmstd_boxplot(subject_grp==2 & dgmstd_grp==3));
    table_pvals(34,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(39).fig = figure(39);
    set(h(39).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dgmstd_grp)/2-0.20 unique(dgmstd_grp)/2+0.20]', repmat(dgmstd_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dgmstd_grp)/2-0.20 unique(dgmstd_grp)/2+0.20]', repmat(dgmstd_median',1,2)','m-','LineWidth',8)
    scatter(dgmstd_grp(subject_grp==3)/2, dgmstd_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dgmstd_grp(subject_grp==1)/2, dgmstd_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dgmstd_grp(subject_grp==4)/2, dgmstd_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dgmstd_grp(subject_grp==2)/2, dgmstd_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.790,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.760,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.670,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.640,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.790,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.790,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.760,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.760,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.670,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.670,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.640,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.640,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d standard deviation value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d standard deviation values','FontSize',18)
    axis([0.5/2 3.5/2 0.00 0.800])
    % axis([0.25 1.75 0.00 0.035])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_std_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==1));
    pCP_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==1));
    pRP_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==1));
    pSC_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==1));
    pSR_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==1));
    pSP_ZOOMint = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==1),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==2),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==2));
    pCR_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==3));
    pCP_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==3));
    pRP_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==3));
    pSC_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==1 & FAWMmean_grp==3));
    pSR_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==3 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==3));
    pSP_RESOLVE = ranksum(FAWMmean_boxplot(subject_grp==4 & FAWMmean_grp==3),FAWMmean_boxplot(subject_grp==2 & FAWMmean_grp==3));
    table_pvals(35,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(40).fig = figure(40);
    set(h(40).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMmean_grp)/2-0.20 unique(FAWMmean_grp)/2+0.20]', repmat(FAWMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMmean_grp)/2-0.20 unique(FAWMmean_grp)/2+0.20]', repmat(FAWMmean_median',1,2)','m-','LineWidth',8)
    scatter(FAWMmean_grp(subject_grp==3)/2, FAWMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMmean_grp(subject_grp==1)/2, FAWMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMmean_grp(subject_grp==4)/2, FAWMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMmean_grp(subject_grp==2)/2, FAWMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.739,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.721,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.712,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.703,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.694,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA mean value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.80])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mean_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==1));
    pCP_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==1));
    pRP_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==1));
    pSC_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==1));
    pSR_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==1));
    pSP_ZOOMint = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==1),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==1));
    pCR_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==2));
    pCP_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==2));
    pRP_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==2));
    pSC_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==2));
    pSR_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==2));
    pSP_ZOOMnotint = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==2),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==2));
    pCR_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==3));
    pCP_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==3));
    pRP_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==3));
    pSC_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==1 & FAGMmean_grp==3));
    pSR_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==3 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==3));
    pSP_RESOLVE = ranksum(FAGMmean_boxplot(subject_grp==4 & FAGMmean_grp==3),FAGMmean_boxplot(subject_grp==2 & FAGMmean_grp==3));
    table_pvals(36,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(41).fig = figure(41);
    set(h(41).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAGMmean_grp)/2-0.20 unique(FAGMmean_grp)/2+0.20]', repmat(FAGMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAGMmean_grp)/2-0.20 unique(FAGMmean_grp)/2+0.20]', repmat(FAGMmean_median',1,2)','m-','LineWidth',8)
    scatter(FAGMmean_grp(subject_grp==3)/2, FAGMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMmean_grp(subject_grp==1)/2, FAGMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMmean_grp(subject_grp==4)/2, FAGMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAGMmean_grp(subject_grp==2)/2, FAGMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.739,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.721,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.712,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.703,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.694,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA mean value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.80])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mean_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MODE VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==1));
    pCP_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==1));
    pRP_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==1));
    pSC_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==1));
    pSR_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==1));
    pSP_ZOOMint = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==1),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==2),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==2));
    pCR_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==3));
    pCP_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==3));
    pRP_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==3));
    pSC_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==1 & FAWMmode_grp==3));
    pSR_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==3 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==3));
    pSP_RESOLVE = ranksum(FAWMmode_boxplot(subject_grp==4 & FAWMmode_grp==3),FAWMmode_boxplot(subject_grp==2 & FAWMmode_grp==3));
    table_pvals(37,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(42).fig = figure(42);
    set(h(42).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMmode_grp)/2-0.20 unique(FAWMmode_grp)/2+0.20]', repmat(FAWMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMmode_grp)/2-0.20 unique(FAWMmode_grp)/2+0.20]', repmat(FAWMmode_median',1,2)','m-','LineWidth',8)
    scatter(FAWMmode_grp(subject_grp==3)/2, FAWMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMmode_grp(subject_grp==1)/2, FAWMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMmode_grp(subject_grp==4)/2, FAWMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMmode_grp(subject_grp==2)/2, FAWMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.940,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.940,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.920,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.920,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.900,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.880,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.880,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.860,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.840,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.840,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.940,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.940,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.920,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.920,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.900,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.880,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.880,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.860,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.840,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.840,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.940,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.940,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.920,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.920,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.900,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.880,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.880,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.860,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.860,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.840,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.840,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA mode value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.95])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mode_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MODE VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==1));
    pCP_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==1));
    pRP_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==1));
    pSC_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==1));
    pSR_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==1));
    pSP_ZOOMint = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==1),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==1));
    pCR_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==2));
    pCP_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==2));
    pRP_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==2));
    pSC_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==2));
    pSR_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==2));
    pSP_ZOOMnotint = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==2),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==2));
    pCR_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==3));
    pCP_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==3));
    pRP_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==3));
    pSC_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==1 & FAGMmode_grp==3));
    pSR_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==3 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==3));
    pSP_RESOLVE = ranksum(FAGMmode_boxplot(subject_grp==4 & FAGMmode_grp==3),FAGMmode_boxplot(subject_grp==2 & FAGMmode_grp==3));
    table_pvals(38,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(43).fig = figure(43);
    set(h(43).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAGMmode_grp)/2-0.20 unique(FAGMmode_grp)/2+0.20]', repmat(FAGMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAGMmode_grp)/2-0.20 unique(FAGMmode_grp)/2+0.20]', repmat(FAGMmode_median',1,2)','m-','LineWidth',8)
    scatter(FAGMmode_grp(subject_grp==3)/2, FAGMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMmode_grp(subject_grp==1)/2, FAGMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMmode_grp(subject_grp==4)/2, FAGMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAGMmode_grp(subject_grp==2)/2, FAGMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.739,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.730,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.721,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.712,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.703,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.694,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.739,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.739,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.730,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.730,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.721,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.721,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.712,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.712,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.703,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.703,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.694,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.694,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA mode value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.4 0.80])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mode_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==1));
    pCP_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==1));
    pRP_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==1));
    pSC_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==1));
    pSR_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==1));
    pSP_ZOOMint = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==1),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==2),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==2));
    pCR_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==3));
    pCP_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==3));
    pRP_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==3));
    pSC_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==1 & f1WMmean_grp==3));
    pSR_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==3 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==3));
    pSP_RESOLVE = ranksum(f1WMmean_boxplot(subject_grp==4 & f1WMmean_grp==3),f1WMmean_boxplot(subject_grp==2 & f1WMmean_grp==3));
    table_pvals(39,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(44).fig = figure(44);
    set(h(44).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMmean_grp)/2-0.20 unique(f1WMmean_grp)/2+0.20]', repmat(f1WMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMmean_grp)/2-0.20 unique(f1WMmean_grp)/2+0.20]', repmat(f1WMmean_median',1,2)','m-','LineWidth',8)
    scatter(f1WMmean_grp(subject_grp==3)/2, f1WMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMmean_grp(subject_grp==1)/2, f1WMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMmean_grp(subject_grp==4)/2, f1WMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMmean_grp(subject_grp==2)/2, f1WMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.630,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.570,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.570,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.540,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.540,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.610,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.610,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.480,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.480,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.630,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.570,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.570,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.540,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.540,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.610,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.610,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.480,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.480,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.630,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.600,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.570,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.570,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.540,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.540,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.610,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.610,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.480,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.480,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 mean value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.25 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mean_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==1));
    pCP_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==1));
    pRP_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==1));
    pSC_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==1));
    pSR_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==1));
    pSP_ZOOMint = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==1),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==1));
    pCR_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==2));
    pCP_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==2));
    pRP_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==2));
    pSC_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==2));
    pSR_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==2));
    pSP_ZOOMnotint = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==2),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==2));
    pCR_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==3));
    pCP_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==3));
    pRP_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==3));
    pSC_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==1 & f1GMmean_grp==3));
    pSR_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==3 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==3));
    pSP_RESOLVE = ranksum(f1GMmean_boxplot(subject_grp==4 & f1GMmean_grp==3),f1GMmean_boxplot(subject_grp==2 & f1GMmean_grp==3));
    table_pvals(40,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(45).fig = figure(45);
    set(h(45).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1GMmean_grp)/2-0.20 unique(f1GMmean_grp)/2+0.20]', repmat(f1GMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1GMmean_grp)/2-0.20 unique(f1GMmean_grp)/2+0.20]', repmat(f1GMmean_median',1,2)','m-','LineWidth',8)
    scatter(f1GMmean_grp(subject_grp==3)/2, f1GMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMmean_grp(subject_grp==1)/2, f1GMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMmean_grp(subject_grp==4)/2, f1GMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1GMmean_grp(subject_grp==2)/2, f1GMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.680,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.655,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.630,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.605,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.580,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.555,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 mean value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.25 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mean_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MODE VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==1));
    pCP_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==1));
    pRP_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==1));
    pSC_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==1));
    pSR_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==1));
    pSP_ZOOMint = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==1),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==2),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==2));
    pCR_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==3));
    pCP_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==3));
    pRP_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==3));
    pSC_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==1 & f1WMmode_grp==3));
    pSR_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==3 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==3));
    pSP_RESOLVE = ranksum(f1WMmode_boxplot(subject_grp==4 & f1WMmode_grp==3),f1WMmode_boxplot(subject_grp==2 & f1WMmode_grp==3));
    table_pvals(41,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(46).fig = figure(46);
    set(h(46).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMmode_grp)/2-0.20 unique(f1WMmode_grp)/2+0.20]', repmat(f1WMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMmode_grp)/2-0.20 unique(f1WMmode_grp)/2+0.20]', repmat(f1WMmode_median',1,2)','m-','LineWidth',8)
    scatter(f1WMmode_grp(subject_grp==3)/2, f1WMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMmode_grp(subject_grp==1)/2, f1WMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMmode_grp(subject_grp==4)/2, f1WMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMmode_grp(subject_grp==2)/2, f1WMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.680,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.655,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.630,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.605,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.580,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.555,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 mode value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.20 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mode_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MODE VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==1));
    pCP_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==1));
    pRP_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==1));
    pSC_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==1));
    pSR_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==1));
    pSP_ZOOMint = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==1),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==1));
    pCR_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==2));
    pCP_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==2));
    pRP_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==2));
    pSC_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==2));
    pSR_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==2));
    pSP_ZOOMnotint = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==2),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==2));
    pCR_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==3));
    pCP_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==3));
    pRP_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==3));
    pSC_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==1 & f1GMmode_grp==3));
    pSR_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==3 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==3));
    pSP_RESOLVE = ranksum(f1GMmode_boxplot(subject_grp==4 & f1GMmode_grp==3),f1GMmode_boxplot(subject_grp==2 & f1GMmode_grp==3));
    table_pvals(42,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(47).fig = figure(47);
    set(h(47).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1GMmode_grp)/2-0.20 unique(f1GMmode_grp)/2+0.20]', repmat(f1GMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1GMmode_grp)/2-0.20 unique(f1GMmode_grp)/2+0.20]', repmat(f1GMmode_median',1,2)','m-','LineWidth',8)
    scatter(f1GMmode_grp(subject_grp==3)/2, f1GMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMmode_grp(subject_grp==1)/2, f1GMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMmode_grp(subject_grp==4)/2, f1GMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1GMmode_grp(subject_grp==2)/2, f1GMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.680,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.655,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.630,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.605,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.580,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.555,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.680,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.680,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.655,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.655,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.630,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.630,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.605,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.605,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.580,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.580,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.555,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.555,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 mode value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.25 0.70])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mode_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MEAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==1));
    pCP_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==1));
    pRP_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==1));
    pSC_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==1));
    pSR_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==1));
    pSP_ZOOMint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==1),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==2),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==2));
    pCR_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==3));
    pCP_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==3));
    pRP_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==3));
    pSC_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==1 & FAWMGMdiffmean_grp==3));
    pSR_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==3 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==3));
    pSP_RESOLVE = ranksum(FAWMGMdiffmean_boxplot(subject_grp==4 & FAWMGMdiffmean_grp==3),FAWMGMdiffmean_boxplot(subject_grp==2 & FAWMGMdiffmean_grp==3));
    table_pvals(43,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(48).fig = figure(48);
    set(h(48).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMGMdiffmean_grp)/2-0.20 unique(FAWMGMdiffmean_grp)/2+0.20]', repmat(FAWMGMdiffmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMGMdiffmean_grp)/2-0.20 unique(FAWMGMdiffmean_grp)/2+0.20]', repmat(FAWMGMdiffmean_median',1,2)','m-','LineWidth',8)
    scatter(FAWMGMdiffmean_grp(subject_grp==3)/2, FAWMGMdiffmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmean_grp(subject_grp==1)/2, FAWMGMdiffmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmean_grp(subject_grp==4)/2, FAWMGMdiffmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmean_grp(subject_grp==2)/2, FAWMGMdiffmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.18,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.16,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.14,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.12,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.08,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mean FA values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA mean differences','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mean_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MEAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==1));
    pCP_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==1));
    pRP_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==1));
    pSC_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==1));
    pSR_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==1));
    pSP_ZOOMint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==1),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==2),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==2));
    pCR_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==3));
    pCP_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==3));
    pRP_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==3));
    pSC_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==1 & f1WMGMdiffmean_grp==3));
    pSR_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==3 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==3));
    pSP_RESOLVE = ranksum(f1WMGMdiffmean_boxplot(subject_grp==4 & f1WMGMdiffmean_grp==3),f1WMGMdiffmean_boxplot(subject_grp==2 & f1WMGMdiffmean_grp==3));
    table_pvals(44,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(49).fig = figure(49);
    set(h(49).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMGMdiffmean_grp)/2-0.20 unique(f1WMGMdiffmean_grp)/2+0.20]', repmat(f1WMGMdiffmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMGMdiffmean_grp)/2-0.20 unique(f1WMGMdiffmean_grp)/2+0.20]', repmat(f1WMGMdiffmean_median',1,2)','m-','LineWidth',8)
    scatter(f1WMGMdiffmean_grp(subject_grp==3)/2, f1WMGMdiffmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmean_grp(subject_grp==1)/2, f1WMGMdiffmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmean_grp(subject_grp==4)/2, f1WMGMdiffmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmean_grp(subject_grp==2)/2, f1WMGMdiffmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.18,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.16,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.14,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.12,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.08,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mean f_1 values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject f_1 mean values','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mean_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==1));
    pCP_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==1));
    pRP_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==1));
    pSC_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==1));
    pSR_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==1));
    pSP_ZOOMint = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==1),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==2),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==2));
    pCR_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==3));
    pCP_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==3));
    pRP_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==3));
    pSC_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==1 & MDWMmean_grp==3));
    pSR_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==3 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==3));
    pSP_RESOLVE = ranksum(MDWMmean_boxplot(subject_grp==4 & MDWMmean_grp==3),MDWMmean_boxplot(subject_grp==2 & MDWMmean_grp==3));
    table_pvals(45,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(50).fig = figure(50);
    set(h(50).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMmean_grp)/2-0.20 unique(MDWMmean_grp)/2+0.20]', repmat(MDWMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMmean_grp)/2-0.20 unique(MDWMmean_grp)/2+0.20]', repmat(MDWMmean_median',1,2)','m-','LineWidth',8)
    scatter(MDWMmean_grp(subject_grp==3)/2, MDWMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMmean_grp(subject_grp==1)/2, MDWMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMmean_grp(subject_grp==4)/2, MDWMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMmean_grp(subject_grp==2)/2, MDWMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD mean value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.7 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mean_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==1));
    pCP_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==1));
    pRP_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==1));
    pSC_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==1));
    pSR_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==1));
    pSP_ZOOMint = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==1),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==1));
    pCR_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==2));
    pCP_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==2));
    pRP_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==2));
    pSC_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==2));
    pSR_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==2));
    pSP_ZOOMnotint = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==2),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==2));
    pCR_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==3));
    pCP_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==3));
    pRP_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==3));
    pSC_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==1 & MDGMmean_grp==3));
    pSR_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==3 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==3));
    pSP_RESOLVE = ranksum(MDGMmean_boxplot(subject_grp==4 & MDGMmean_grp==3),MDGMmean_boxplot(subject_grp==2 & MDGMmean_grp==3));
    table_pvals(46,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(51).fig = figure(51);
    set(h(51).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDGMmean_grp)/2-0.20 unique(MDGMmean_grp)/2+0.20]', repmat(MDGMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDGMmean_grp)/2-0.20 unique(MDGMmean_grp)/2+0.20]', repmat(MDGMmean_median',1,2)','m-','LineWidth',8)
    scatter(MDGMmean_grp(subject_grp==3)/2, MDGMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMmean_grp(subject_grp==1)/2, MDGMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMmean_grp(subject_grp==4)/2, MDGMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDGMmean_grp(subject_grp==2)/2, MDGMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.570,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.570,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.520,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.520,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.470,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.470,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.420,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.420,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD mean value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD mean values','FontSize',18)
    axis([0.5/2 3.5/2 0.7 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mean_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MODE VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==1));
    pCP_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==1));
    pRP_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==1));
    pSC_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==1));
    pSR_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==1));
    pSP_ZOOMint = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==1),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==2),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==2));
    pCR_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==3));
    pCP_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==3));
    pRP_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==3));
    pSC_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==1 & MDWMmode_grp==3));
    pSR_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==3 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==3));
    pSP_RESOLVE = ranksum(MDWMmode_boxplot(subject_grp==4 & MDWMmode_grp==3),MDWMmode_boxplot(subject_grp==2 & MDWMmode_grp==3));
    table_pvals(47,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(52).fig = figure(52);
    set(h(52).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMmode_grp)/2-0.20 unique(MDWMmode_grp)/2+0.20]', repmat(MDWMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMmode_grp)/2-0.20 unique(MDWMmode_grp)/2+0.20]', repmat(MDWMmode_median',1,2)','m-','LineWidth',8)
    scatter(MDWMmode_grp(subject_grp==3)/2, MDWMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMmode_grp(subject_grp==1)/2, MDWMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMmode_grp(subject_grp==4)/2, MDWMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMmode_grp(subject_grp==2)/2, MDWMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD mode value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.7 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mode_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MODE VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==1));
    pCP_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==1));
    pRP_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==1));
    pSC_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==1));
    pSR_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==1));
    pSP_ZOOMint = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==1),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==1));
    pCR_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==2));
    pCP_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==2));
    pRP_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==2));
    pSC_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==2));
    pSR_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==2));
    pSP_ZOOMnotint = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==2),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==2));
    pCR_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==3));
    pCP_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==3));
    pRP_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==3));
    pSC_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==1 & MDGMmode_grp==3));
    pSR_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==3 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==3));
    pSP_RESOLVE = ranksum(MDGMmode_boxplot(subject_grp==4 & MDGMmode_grp==3),MDGMmode_boxplot(subject_grp==2 & MDGMmode_grp==3));
    table_pvals(48,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(53).fig = figure(53);
    set(h(53).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDGMmode_grp)/2-0.20 unique(MDGMmode_grp)/2+0.20]', repmat(MDGMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDGMmode_grp)/2-0.20 unique(MDGMmode_grp)/2+0.20]', repmat(MDGMmode_median',1,2)','m-','LineWidth',8)
    scatter(MDGMmode_grp(subject_grp==3)/2, MDGMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMmode_grp(subject_grp==1)/2, MDGMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMmode_grp(subject_grp==4)/2, MDGMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDGMmode_grp(subject_grp==2)/2, MDGMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,1.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,1.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,1.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,1.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,1.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,1.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,1.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,1.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,1.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,1.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD mode value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.6 1.6])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mode_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MEAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==1));
    pCP_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==1));
    pRP_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==1));
    pSC_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==1));
    pSR_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==1));
    pSP_ZOOMint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==1),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==2),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==2));
    pCR_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==3));
    pCP_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==3));
    pRP_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==3));
    pSC_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==1 & MDWMGMdiffmean_grp==3));
    pSR_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==3 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==3));
    pSP_RESOLVE = ranksum(MDWMGMdiffmean_boxplot(subject_grp==4 & MDWMGMdiffmean_grp==3),MDWMGMdiffmean_boxplot(subject_grp==2 & MDWMGMdiffmean_grp==3));
    table_pvals(49,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(54).fig = figure(54);
    set(h(54).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMGMdiffmean_grp)/2-0.20 unique(MDWMGMdiffmean_grp)/2+0.20]', repmat(MDWMGMdiffmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMGMdiffmean_grp)/2-0.20 unique(MDWMGMdiffmean_grp)/2+0.20]', repmat(MDWMGMdiffmean_median',1,2)','m-','LineWidth',8)
    scatter(MDWMGMdiffmean_grp(subject_grp==3)/2, MDWMGMdiffmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmean_grp(subject_grp==1)/2, MDWMGMdiffmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmean_grp(subject_grp==4)/2, MDWMGMdiffmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmean_grp(subject_grp==2)/2, MDWMGMdiffmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.49,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.47,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.45,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.43,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.41,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.39,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mean MD values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD mean diferences','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.50])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mean_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEAN VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==3 & dWMmean_grp==1));
    pCP_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==1));
    pRP_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==1));
    pSC_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==1 & dWMmean_grp==1));
    pSR_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==4 & dWMmean_grp==1));
    pSP_ZOOMint = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==1),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==1));
    pCR_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==3 & dWMmean_grp==2));
    pCP_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==2));
    pRP_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==2));
    pSC_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==1 & dWMmean_grp==2));
    pSR_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==4 & dWMmean_grp==2));
    pSP_ZOOMnotint = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==2),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==2));
    pCR_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==3 & dWMmean_grp==3));
    pCP_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==1 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==3));
    pRP_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==3));
    pSC_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==1 & dWMmean_grp==3));
    pSR_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==3 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==4 & dWMmean_grp==3));
    pSP_RESOLVE = ranksum(dWMmean_boxplot(subject_grp==4 & dWMmean_grp==3),dWMmean_boxplot(subject_grp==2 & dWMmean_grp==3));
    table_pvals(50,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(55).fig = figure(55);
    set(h(55).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMmean_grp)/2-0.20 unique(dWMmean_grp)/2+0.20]', repmat(dWMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMmean_grp)/2-0.20 unique(dWMmean_grp)/2+0.20]', repmat(dWMmean_median',1,2)','m-','LineWidth',8)
    scatter(dWMmean_grp(subject_grp==3)/2, dWMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMmean_grp(subject_grp==1)/2, dWMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMmean_grp(subject_grp==4)/2, dWMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMmean_grp(subject_grp==2)/2, dWMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d mean value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d mean values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mean_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEAN VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==3 & dGMmean_grp==1));
    pCP_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==1));
    pRP_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==1));
    pSC_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==1 & dGMmean_grp==1));
    pSR_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==4 & dGMmean_grp==1));
    pSP_ZOOMint = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==1),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==1));
    pCR_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==3 & dGMmean_grp==2));
    pCP_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==2));
    pRP_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==2));
    pSC_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==1 & dGMmean_grp==2));
    pSR_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==4 & dGMmean_grp==2));
    pSP_ZOOMnotint = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==2),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==2));
    pCR_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==3 & dGMmean_grp==3));
    pCP_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==1 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==3));
    pRP_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==3));
    pSC_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==1 & dGMmean_grp==3));
    pSR_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==3 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==4 & dGMmean_grp==3));
    pSP_RESOLVE = ranksum(dGMmean_boxplot(subject_grp==4 & dGMmean_grp==3),dGMmean_boxplot(subject_grp==2 & dGMmean_grp==3));
    table_pvals(51,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(56).fig = figure(56);
    set(h(56).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dGMmean_grp)/2-0.20 unique(dGMmean_grp)/2+0.20]', repmat(dGMmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dGMmean_grp)/2-0.20 unique(dGMmean_grp)/2+0.20]', repmat(dGMmean_median',1,2)','m-','LineWidth',8)
    scatter(dGMmean_grp(subject_grp==3)/2, dGMmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMmean_grp(subject_grp==1)/2, dGMmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMmean_grp(subject_grp==4)/2, dGMmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dGMmean_grp(subject_grp==2)/2, dGMmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.510,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.510,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.390,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.390,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d mean value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d mean values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mean_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MEAN WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==1));
    pCP_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==1));
    pRP_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==1));
    pSC_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==1));
    pSR_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==1));
    pSP_ZOOMint = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==1),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==1));
    pCR_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==2));
    pCP_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==2));
    pRP_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==2));
    pSC_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==2));
    pSR_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==2));
    pSP_ZOOMnotint = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==2),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==2));
    pCR_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==3));
    pCP_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==3));
    pRP_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==3));
    pSC_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==1 & dWMGMdiffmean_grp==3));
    pSR_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==3 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==3));
    pSP_RESOLVE = ranksum(dWMGMdiffmean_boxplot(subject_grp==4 & dWMGMdiffmean_grp==3),dWMGMdiffmean_boxplot(subject_grp==2 & dWMGMdiffmean_grp==3));
    table_pvals(52,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(57).fig = figure(57);
    set(h(57).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMGMdiffmean_grp)/2-0.20 unique(dWMGMdiffmean_grp)/2+0.20]', repmat(dWMGMdiffmean_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMGMdiffmean_grp)/2-0.20 unique(dWMGMdiffmean_grp)/2+0.20]', repmat(dWMGMdiffmean_median',1,2)','m-','LineWidth',8)
    scatter(dWMGMdiffmean_grp(subject_grp==3)/2, dWMGMdiffmean_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmean_grp(subject_grp==1)/2, dWMGMdiffmean_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmean_grp(subject_grp==4)/2, dWMGMdiffmean_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmean_grp(subject_grp==2)/2, dWMGMdiffmean_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.83,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.80,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.77,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.74,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.71,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.68,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mean d values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d mean differences','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.85])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mean_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MODE VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==3 & dWMmode_grp==1));
    pCP_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==1));
    pRP_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==1));
    pSC_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==1 & dWMmode_grp==1));
    pSR_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==4 & dWMmode_grp==1));
    pSP_ZOOMint = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==1),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==1));
    pCR_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==3 & dWMmode_grp==2));
    pCP_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==2));
    pRP_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==2));
    pSC_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==1 & dWMmode_grp==2));
    pSR_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==4 & dWMmode_grp==2));
    pSP_ZOOMnotint = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==2),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==2));
    pCR_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==3 & dWMmode_grp==3));
    pCP_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==1 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==3));
    pRP_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==3));
    pSC_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==1 & dWMmode_grp==3));
    pSR_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==3 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==4 & dWMmode_grp==3));
    pSP_RESOLVE = ranksum(dWMmode_boxplot(subject_grp==4 & dWMmode_grp==3),dWMmode_boxplot(subject_grp==2 & dWMmode_grp==3));
    table_pvals(53,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(58).fig = figure(58);
    set(h(58).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMmode_grp)/2-0.20 unique(dWMmode_grp)/2+0.20]', repmat(dWMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMmode_grp)/2-0.20 unique(dWMmode_grp)/2+0.20]', repmat(dWMmode_median',1,2)','m-','LineWidth',8)
    scatter(dWMmode_grp(subject_grp==3)/2, dWMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMmode_grp(subject_grp==1)/2, dWMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMmode_grp(subject_grp==4)/2, dWMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMmode_grp(subject_grp==2)/2, dWMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d mode value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d mode values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mode_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MODE VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==3 & dGMmode_grp==1));
    pCP_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==1));
    pRP_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==1));
    pSC_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==1 & dGMmode_grp==1));
    pSR_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==4 & dGMmode_grp==1));
    pSP_ZOOMint = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==1),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==1));
    pCR_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==3 & dGMmode_grp==2));
    pCP_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==2));
    pRP_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==2));
    pSC_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==1 & dGMmode_grp==2));
    pSR_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==4 & dGMmode_grp==2));
    pSP_ZOOMnotint = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==2),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==2));
    pCR_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==3 & dGMmode_grp==3));
    pCP_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==1 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==3));
    pRP_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==3));
    pSC_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==1 & dGMmode_grp==3));
    pSR_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==3 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==4 & dGMmode_grp==3));
    pSP_RESOLVE = ranksum(dGMmode_boxplot(subject_grp==4 & dGMmode_grp==3),dGMmode_boxplot(subject_grp==2 & dGMmode_grp==3));
    table_pvals(54,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(59).fig = figure(59);
    set(h(59).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dGMmode_grp)/2-0.20 unique(dGMmode_grp)/2+0.20]', repmat(dGMmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dGMmode_grp)/2-0.20 unique(dGMmode_grp)/2+0.20]', repmat(dGMmode_median',1,2)','m-','LineWidth',8)
    scatter(dGMmode_grp(subject_grp==3)/2, dGMmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMmode_grp(subject_grp==1)/2, dGMmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMmode_grp(subject_grp==4)/2, dGMmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dGMmode_grp(subject_grp==2)/2, dGMmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,2.550,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.450,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.400,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.350,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.300,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,2.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,2.550,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.450,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.400,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.350,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.350,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.300,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d mode value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.9 2.65])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mode_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA SKEWNESS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==1));
    pCP_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==1));
    pRP_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==1));
    pSC_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==1));
    pSR_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==1));
    pSP_ZOOMint = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==1),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==2),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==2));
    pCR_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==3));
    pCP_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==3));
    pRP_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==3));
    pSC_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==1 & FAWMskewness_grp==3));
    pSR_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==3 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==3));
    pSP_RESOLVE = ranksum(FAWMskewness_boxplot(subject_grp==4 & FAWMskewness_grp==3),FAWMskewness_boxplot(subject_grp==2 & FAWMskewness_grp==3));
    table_pvals(55,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(60).fig = figure(60);
    set(h(60).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMskewness_grp)/2-0.20 unique(FAWMskewness_grp)/2+0.20]', repmat(FAWMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMskewness_grp)/2-0.20 unique(FAWMskewness_grp)/2+0.20]', repmat(FAWMskewness_median',1,2)','m-','LineWidth',8)
    scatter(FAWMskewness_grp(subject_grp==3)/2, FAWMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMskewness_grp(subject_grp==1)/2, FAWMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMskewness_grp(subject_grp==4)/2, FAWMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMskewness_grp(subject_grp==2)/2, FAWMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.650,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.600,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.550,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.500,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.450,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.400,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA skewness value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -0.9 0.7])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_skewness_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA SKEWNESS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==1));
    pCP_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==1));
    pRP_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==1));
    pSC_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==1));
    pSR_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==1));
    pSP_ZOOMint = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==1),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==2),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==2));
    pCR_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==3));
    pCP_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==3));
    pRP_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==3));
    pSC_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==1 & FAGMskewness_grp==3));
    pSR_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==3 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==3));
    pSP_RESOLVE = ranksum(FAGMskewness_boxplot(subject_grp==4 & FAGMskewness_grp==3),FAGMskewness_boxplot(subject_grp==2 & FAGMskewness_grp==3));
    table_pvals(56,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(61).fig = figure(61);
    set(h(61).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAGMskewness_grp)/2-0.20 unique(FAGMskewness_grp)/2+0.20]', repmat(FAGMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAGMskewness_grp)/2-0.20 unique(FAGMskewness_grp)/2+0.20]', repmat(FAGMskewness_median',1,2)','m-','LineWidth',8)
    scatter(FAGMskewness_grp(subject_grp==3)/2, FAGMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMskewness_grp(subject_grp==1)/2, FAGMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMskewness_grp(subject_grp==4)/2, FAGMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAGMskewness_grp(subject_grp==2)/2, FAGMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.650,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.600,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.550,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.500,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.450,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.400,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA skewness value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -0.9 0.7])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_skewness_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 SKEWNESS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==1));
    pCP_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==1));
    pRP_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==1));
    pSC_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==1));
    pSR_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==1));
    pSP_ZOOMint = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==1),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==2),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==2));
    pCR_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==3));
    pCP_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==3));
    pRP_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==3));
    pSC_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==1 & f1WMskewness_grp==3));
    pSR_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==3 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==3));
    pSP_RESOLVE = ranksum(f1WMskewness_boxplot(subject_grp==4 & f1WMskewness_grp==3),f1WMskewness_boxplot(subject_grp==2 & f1WMskewness_grp==3));
    table_pvals(57,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(62).fig = figure(62);
    set(h(62).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMskewness_grp)/2-0.20 unique(f1WMskewness_grp)/2+0.20]', repmat(f1WMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMskewness_grp)/2-0.20 unique(f1WMskewness_grp)/2+0.20]', repmat(f1WMskewness_median',1,2)','m-','LineWidth',8)
    scatter(f1WMskewness_grp(subject_grp==3)/2, f1WMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMskewness_grp(subject_grp==1)/2, f1WMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMskewness_grp(subject_grp==4)/2, f1WMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMskewness_grp(subject_grp==2)/2, f1WMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.650,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.600,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.550,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.500,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.450,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.400,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 skewness value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -0.9 0.7])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_skewness_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 SKEWNESS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==1));
    pCP_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==1));
    pRP_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==1));
    pSC_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==1));
    pSR_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==1));
    pSP_ZOOMint = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==1),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==2),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==2));
    pCR_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==3));
    pCP_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==3));
    pRP_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==3));
    pSC_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==1 & f1GMskewness_grp==3));
    pSR_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==3 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==3));
    pSP_RESOLVE = ranksum(f1GMskewness_boxplot(subject_grp==4 & f1GMskewness_grp==3),f1GMskewness_boxplot(subject_grp==2 & f1GMskewness_grp==3));
    table_pvals(58,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(63).fig = figure(63);
    set(h(63).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1GMskewness_grp)/2-0.20 unique(f1GMskewness_grp)/2+0.20]', repmat(f1GMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1GMskewness_grp)/2-0.20 unique(f1GMskewness_grp)/2+0.20]', repmat(f1GMskewness_median',1,2)','m-','LineWidth',8)
    scatter(f1GMskewness_grp(subject_grp==3)/2, f1GMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMskewness_grp(subject_grp==1)/2, f1GMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMskewness_grp(subject_grp==4)/2, f1GMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1GMskewness_grp(subject_grp==2)/2, f1GMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.650,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.600,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.550,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.500,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.450,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.400,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.650,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.650,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.600,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.550,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.550,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.500,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.450,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.450,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.400,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.400,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 skewness value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -0.9 0.7])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_skewness_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD SKEWNESS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==1));
    pCP_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==1));
    pRP_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==1));
    pSC_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==1));
    pSR_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==1));
    pSP_ZOOMint = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==1),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==2),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==2));
    pCR_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==3));
    pCP_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==3));
    pRP_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==3));
    pSC_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==1 & MDWMskewness_grp==3));
    pSR_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==3 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==3));
    pSP_RESOLVE = ranksum(MDWMskewness_boxplot(subject_grp==4 & MDWMskewness_grp==3),MDWMskewness_boxplot(subject_grp==2 & MDWMskewness_grp==3));
    table_pvals(59,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(64).fig = figure(64);
    set(h(64).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMskewness_grp)/2-0.20 unique(MDWMskewness_grp)/2+0.20]', repmat(MDWMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMskewness_grp)/2-0.20 unique(MDWMskewness_grp)/2+0.20]', repmat(MDWMskewness_median',1,2)','m-','LineWidth',8)
    scatter(MDWMskewness_grp(subject_grp==3)/2, MDWMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMskewness_grp(subject_grp==1)/2, MDWMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMskewness_grp(subject_grp==4)/2, MDWMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMskewness_grp(subject_grp==2)/2, MDWMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,3.000,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.900,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.800,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.600,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD skewness value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -1.0 3.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_skewness_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD SKEWNESS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==1));
    pCP_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==1));
    pRP_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==1));
    pSC_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==1));
    pSR_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==1));
    pSP_ZOOMint = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==1),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==2),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==2));
    pCR_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==3));
    pCP_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==3));
    pRP_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==3));
    pSC_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==1 & MDGMskewness_grp==3));
    pSR_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==3 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==3));
    pSP_RESOLVE = ranksum(MDGMskewness_boxplot(subject_grp==4 & MDGMskewness_grp==3),MDGMskewness_boxplot(subject_grp==2 & MDGMskewness_grp==3));
    table_pvals(60,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(65).fig = figure(65);
    set(h(65).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDGMskewness_grp)/2-0.20 unique(MDGMskewness_grp)/2+0.20]', repmat(MDGMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDGMskewness_grp)/2-0.20 unique(MDGMskewness_grp)/2+0.20]', repmat(MDGMskewness_median',1,2)','m-','LineWidth',8)
    scatter(MDGMskewness_grp(subject_grp==3)/2, MDGMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMskewness_grp(subject_grp==1)/2, MDGMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMskewness_grp(subject_grp==4)/2, MDGMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDGMskewness_grp(subject_grp==2)/2, MDGMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.800,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,3.000,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.900,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.800,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.600,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD skewness value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -1.0 3.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_skewness_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d SKEWNESS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==1));
    pCP_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==1));
    pRP_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==1));
    pSC_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==1));
    pSR_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==1));
    pSP_ZOOMint = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==1),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==2),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==2));
    pCR_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==3));
    pCP_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==3));
    pRP_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==3));
    pSC_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==1 & dWMskewness_grp==3));
    pSR_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==3 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==3));
    pSP_RESOLVE = ranksum(dWMskewness_boxplot(subject_grp==4 & dWMskewness_grp==3),dWMskewness_boxplot(subject_grp==2 & dWMskewness_grp==3));
    table_pvals(61,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(66).fig = figure(66);
    set(h(66).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMskewness_grp)/2-0.20 unique(dWMskewness_grp)/2+0.20]', repmat(dWMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMskewness_grp)/2-0.20 unique(dWMskewness_grp)/2+0.20]', repmat(dWMskewness_median',1,2)','m-','LineWidth',8)
    scatter(dWMskewness_grp(subject_grp==3)/2, dWMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMskewness_grp(subject_grp==1)/2, dWMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMskewness_grp(subject_grp==4)/2, dWMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMskewness_grp(subject_grp==2)/2, dWMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,2.900,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,3.000,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.900,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.800,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.600,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d skewness value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -1.0 3.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_skewness_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d SKEWNESS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==1));
    pCP_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==1));
    pRP_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==1));
    pSC_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==1));
    pSR_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==1));
    pSP_ZOOMint = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==1),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==1));
    pCR_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==2));
    pCP_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==2));
    pRP_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==2));
    pSC_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==2));
    pSR_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==2));
    pSP_ZOOMnotint = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==2),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==2));
    pCR_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==3));
    pCP_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==3));
    pRP_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==3));
    pSC_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==1 & dGMskewness_grp==3));
    pSR_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==3 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==3));
    pSP_RESOLVE = ranksum(dGMskewness_boxplot(subject_grp==4 & dGMskewness_grp==3),dGMskewness_boxplot(subject_grp==2 & dGMskewness_grp==3));
    table_pvals(62,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(67).fig = figure(67);
    set(h(67).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dGMskewness_grp)/2-0.20 unique(dGMskewness_grp)/2+0.20]', repmat(dGMskewness_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dGMskewness_grp)/2-0.20 unique(dGMskewness_grp)/2+0.20]', repmat(dGMskewness_median',1,2)','m-','LineWidth',8)
    scatter(dGMskewness_grp(subject_grp==3)/2, dGMskewness_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMskewness_grp(subject_grp==1)/2, dGMskewness_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMskewness_grp(subject_grp==4)/2, dGMskewness_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dGMskewness_grp(subject_grp==2)/2, dGMskewness_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,3.000,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.900,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.800,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.700,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.600,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,2.500,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,3.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,3.000,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,2.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.900,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,2.800,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.800,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,2.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.700,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,2.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.600,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,2.500,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,2.500,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'d skewness value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d skewness values','FontSize',18)
    axis([0.5/2 3.5/2 -1.0 3.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_skewness_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA KURTOSIS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==1),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==2),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==1 & FAWMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==3 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(FAWMkurtosis_boxplot(subject_grp==4 & FAWMkurtosis_grp==3),FAWMkurtosis_boxplot(subject_grp==2 & FAWMkurtosis_grp==3));
    table_pvals(63,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(68).fig = figure(68);
    set(h(68).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMkurtosis_grp)/2-0.20 unique(FAWMkurtosis_grp)/2+0.20]', repmat(FAWMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMkurtosis_grp)/2-0.20 unique(FAWMkurtosis_grp)/2+0.20]', repmat(FAWMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(FAWMkurtosis_grp(subject_grp==3)/2, FAWMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMkurtosis_grp(subject_grp==1)/2, FAWMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMkurtosis_grp(subject_grp==4)/2, FAWMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMkurtosis_grp(subject_grp==2)/2, FAWMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,4.250,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.150,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.050,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.950,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.850,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.750,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA kurtosis value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 2.0 4.3])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_kurtosis_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA KURTOSIS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==1),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==2),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==1 & FAGMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==3 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(FAGMkurtosis_boxplot(subject_grp==4 & FAGMkurtosis_grp==3),FAGMkurtosis_boxplot(subject_grp==2 & FAGMkurtosis_grp==3));
    table_pvals(64,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(69).fig = figure(69);
    set(h(69).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAGMkurtosis_grp)/2-0.20 unique(FAGMkurtosis_grp)/2+0.20]', repmat(FAGMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAGMkurtosis_grp)/2-0.20 unique(FAGMkurtosis_grp)/2+0.20]', repmat(FAGMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(FAGMkurtosis_grp(subject_grp==3)/2, FAGMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMkurtosis_grp(subject_grp==1)/2, FAGMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAGMkurtosis_grp(subject_grp==4)/2, FAGMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAGMkurtosis_grp(subject_grp==2)/2, FAGMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,4.250,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.150,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.050,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.950,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.850,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.750,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'FA kurtosis value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject FA kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 2.0 4.3])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_kurtosis_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 KURTOSIS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==1),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==2),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==1 & f1WMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==3 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(f1WMkurtosis_boxplot(subject_grp==4 & f1WMkurtosis_grp==3),f1WMkurtosis_boxplot(subject_grp==2 & f1WMkurtosis_grp==3));
    table_pvals(65,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(70).fig = figure(70);
    set(h(70).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMkurtosis_grp)/2-0.20 unique(f1WMkurtosis_grp)/2+0.20]', repmat(f1WMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMkurtosis_grp)/2-0.20 unique(f1WMkurtosis_grp)/2+0.20]', repmat(f1WMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(f1WMkurtosis_grp(subject_grp==3)/2, f1WMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMkurtosis_grp(subject_grp==1)/2, f1WMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMkurtosis_grp(subject_grp==4)/2, f1WMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMkurtosis_grp(subject_grp==2)/2, f1WMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,4.250,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.150,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.050,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.950,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.850,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.750,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 kurtosis value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 1.8 4.3])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_kurtosis_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 KURTOSIS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==1),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==2),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==1 & f1GMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==3 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(f1GMkurtosis_boxplot(subject_grp==4 & f1GMkurtosis_grp==3),f1GMkurtosis_boxplot(subject_grp==2 & f1GMkurtosis_grp==3));
    table_pvals(66,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(71).fig = figure(71);
    set(h(71).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1GMkurtosis_grp)/2-0.20 unique(f1GMkurtosis_grp)/2+0.20]', repmat(f1GMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1GMkurtosis_grp)/2-0.20 unique(f1GMkurtosis_grp)/2+0.20]', repmat(f1GMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(f1GMkurtosis_grp(subject_grp==3)/2, f1GMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMkurtosis_grp(subject_grp==1)/2, f1GMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1GMkurtosis_grp(subject_grp==4)/2, f1GMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1GMkurtosis_grp(subject_grp==2)/2, f1GMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,4.250,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.150,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,4.050,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.950,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.850,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,3.750,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,4.250,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,4.250,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,4.150,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.150,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,4.050,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,4.050,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,3.950,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.950,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,3.850,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.850,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,3.750,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,3.750,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'f_1 kurtosis value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject f_1 kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 1.8 4.3])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_kurtosis_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD KURTOSIS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==1),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==2),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==1 & MDWMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==3 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(MDWMkurtosis_boxplot(subject_grp==4 & MDWMkurtosis_grp==3),MDWMkurtosis_boxplot(subject_grp==2 & MDWMkurtosis_grp==3));
    table_pvals(67,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(72).fig = figure(72);
    set(h(72).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMkurtosis_grp)/2-0.20 unique(MDWMkurtosis_grp)/2+0.20]', repmat(MDWMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMkurtosis_grp)/2-0.20 unique(MDWMkurtosis_grp)/2+0.20]', repmat(MDWMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(MDWMkurtosis_grp(subject_grp==3)/2, MDWMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMkurtosis_grp(subject_grp==1)/2, MDWMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMkurtosis_grp(subject_grp==4)/2, MDWMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMkurtosis_grp(subject_grp==2)/2, MDWMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,11.30,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,11.30,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,11.00,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,11.00,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,10.60,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,10.60,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,10.40,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,10.40,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,10.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,10.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,9.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,9.80,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,11.30,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,11.30,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,11.00,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,11.00,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,10.70,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,10.70,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,10.40,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,10.40,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,10.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,10.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,9.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,9.80,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,11.30,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,11.30,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,11.00,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,11.00,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,10.70,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,10.70,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,10.40,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,10.40,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,10.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,10.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,9.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,9.80,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD kurtosis value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 2.0 11.5])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_kurtosis_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD KURTOSIS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==1),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==2),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==1 & MDGMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==3 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(MDGMkurtosis_boxplot(subject_grp==4 & MDGMkurtosis_grp==3),MDGMkurtosis_boxplot(subject_grp==2 & MDGMkurtosis_grp==3));
    table_pvals(68,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(73).fig = figure(73);
    set(h(73).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDGMkurtosis_grp)/2-0.20 unique(MDGMkurtosis_grp)/2+0.20]', repmat(MDGMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDGMkurtosis_grp)/2-0.20 unique(MDGMkurtosis_grp)/2+0.20]', repmat(MDGMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(MDGMkurtosis_grp(subject_grp==3)/2, MDGMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMkurtosis_grp(subject_grp==1)/2, MDGMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDGMkurtosis_grp(subject_grp==4)/2, MDGMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDGMkurtosis_grp(subject_grp==2)/2, MDGMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,17.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,17.700,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,17.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,17.000,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,16.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,16.300,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,15.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,15.600,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,14.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,14.900,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,14.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,14.200,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,17.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,17.700,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,17.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,17.000,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,16.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,16.300,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,15.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,15.600,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,14.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,14.900,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,14.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,14.200,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,17.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,17.700,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,17.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,17.000,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,16.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,16.300,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,15.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,15.600,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,14.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,14.900,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,14.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,14.200,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'MD kurtosis value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject MD kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 2.0 18.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_kurtosis_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d KURTOSIS VALUE DISTRIBUTIONS FROM WM
    pCR_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==1),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==2),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==1 & dWMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==3 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(dWMkurtosis_boxplot(subject_grp==4 & dWMkurtosis_grp==3),dWMkurtosis_boxplot(subject_grp==2 & dWMkurtosis_grp==3));
    table_pvals(69,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(74).fig = figure(74);
    set(h(74).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMkurtosis_grp)/2-0.20 unique(dWMkurtosis_grp)/2+0.20]', repmat(dWMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMkurtosis_grp)/2-0.20 unique(dWMkurtosis_grp)/2+0.20]', repmat(dWMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(dWMkurtosis_grp(subject_grp==3)/2, dWMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMkurtosis_grp(subject_grp==1)/2, dWMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMkurtosis_grp(subject_grp==4)/2, dWMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMkurtosis_grp(subject_grp==2)/2, dWMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,12.700,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,12.000,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,11.300,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,10.600,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,9.900,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,9.200,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,12.700,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,12.000,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,11.300,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,10.600,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,9.900,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,9.200,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,12.700,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,12.000,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,11.300,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,10.600,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,9.900,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,9.200,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end 
    hold off
    title({'d kurtosis value distribution from WM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 13.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_kurtosis_WM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d KURTOSIS VALUE DISTRIBUTIONS FROM GM
    pCR_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==1));
    pCP_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==1));
    pRP_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==1));
    pSC_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==1));
    pSR_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==1));
    pSP_ZOOMint = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==1),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==1));
    pCR_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==2));
    pCP_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==2));
    pRP_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==2));
    pSC_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==2));
    pSR_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==2));
    pSP_ZOOMnotint = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==2),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==2));
    pCR_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==3));
    pCP_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==3));
    pRP_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==3));
    pSC_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==1 & dGMkurtosis_grp==3));
    pSR_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==3 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==3));
    pSP_RESOLVE = ranksum(dGMkurtosis_boxplot(subject_grp==4 & dGMkurtosis_grp==3),dGMkurtosis_boxplot(subject_grp==2 & dGMkurtosis_grp==3));
    table_pvals(70,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(75).fig = figure(75);
    set(h(75).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dGMkurtosis_grp)/2-0.20 unique(dGMkurtosis_grp)/2+0.20]', repmat(dGMkurtosis_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dGMkurtosis_grp)/2-0.20 unique(dGMkurtosis_grp)/2+0.20]', repmat(dGMkurtosis_median',1,2)','m-','LineWidth',8)
    scatter(dGMkurtosis_grp(subject_grp==3)/2, dGMkurtosis_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMkurtosis_grp(subject_grp==1)/2, dGMkurtosis_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dGMkurtosis_grp(subject_grp==4)/2, dGMkurtosis_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dGMkurtosis_grp(subject_grp==2)/2, dGMkurtosis_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,12.700,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,12.000,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,11.300,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,10.600,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,9.900,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,9.200,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,12.700,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,12.000,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,11.300,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,10.600,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,9.900,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,9.200,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,12.700,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,12.700,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,12.000,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,12.000,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,11.300,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,11.300,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,10.600,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,10.600,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,9.900,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,9.900,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,9.200,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,9.200,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end 
    hold off
    title({'d kurtosis value distribution from GM';'C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Single-subject d kurtosis values','FontSize',18)
    axis([0.5/2 3.5/2 1.0 13.0])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_kurtosis_GM_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF FA MODE WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==1));
    pCP_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==1));
    pRP_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==1));
    pSC_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==1));
    pSR_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==1));
    pSP_ZOOMint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==1),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==1));
    pCR_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==2));
    pCP_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==2));
    pRP_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==2));
    pSC_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==2));
    pSR_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==2));
    pSP_ZOOMnotint = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==2),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==2));
    pCR_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==3));
    pCP_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==3));
    pRP_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==3));
    pSC_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==1 & FAWMGMdiffmode_grp==3));
    pSR_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==3 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==3));
    pSP_RESOLVE = ranksum(FAWMGMdiffmode_boxplot(subject_grp==4 & FAWMGMdiffmode_grp==3),FAWMGMdiffmode_boxplot(subject_grp==2 & FAWMGMdiffmode_grp==3));
    table_pvals(71,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(76).fig = figure(76);
    set(h(76).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(FAWMGMdiffmode_grp)/2-0.20 unique(FAWMGMdiffmode_grp)/2+0.20]', repmat(FAWMGMdiffmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(FAWMGMdiffmode_grp)/2-0.20 unique(FAWMGMdiffmode_grp)/2+0.20]', repmat(FAWMGMdiffmode_median',1,2)','m-','LineWidth',8)
    scatter(FAWMGMdiffmode_grp(subject_grp==3)/2, FAWMGMdiffmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmode_grp(subject_grp==1)/2, FAWMGMdiffmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmode_grp(subject_grp==4)/2, FAWMGMdiffmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(FAWMGMdiffmode_grp(subject_grp==2)/2, FAWMGMdiffmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.18,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.16,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.14,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.12,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.08,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mode FA values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject FA mode values','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_FA_mode_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF f1 MODE WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==1));
    pCP_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==1));
    pRP_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==1));
    pSC_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==1));
    pSR_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==1));
    pSP_ZOOMint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==1),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==1));
    pCR_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==2));
    pCP_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==2));
    pRP_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==2));
    pSC_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==2));
    pSR_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==2));
    pSP_ZOOMnotint = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==2),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==2));
    pCR_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==3));
    pCP_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==3));
    pRP_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==3));
    pSC_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==1 & f1WMGMdiffmode_grp==3));
    pSR_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==3 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==3));
    pSP_RESOLVE = ranksum(f1WMGMdiffmode_boxplot(subject_grp==4 & f1WMGMdiffmode_grp==3),f1WMGMdiffmode_boxplot(subject_grp==2 & f1WMGMdiffmode_grp==3));
    table_pvals(72,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(77).fig = figure(77);
    set(h(77).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(f1WMGMdiffmode_grp)/2-0.20 unique(f1WMGMdiffmode_grp)/2+0.20]', repmat(f1WMGMdiffmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(f1WMGMdiffmode_grp)/2-0.20 unique(f1WMGMdiffmode_grp)/2+0.20]', repmat(f1WMGMdiffmode_median',1,2)','m-','LineWidth',8)
    scatter(f1WMGMdiffmode_grp(subject_grp==3)/2, f1WMGMdiffmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmode_grp(subject_grp==1)/2, f1WMGMdiffmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmode_grp(subject_grp==4)/2, f1WMGMdiffmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(f1WMGMdiffmode_grp(subject_grp==2)/2, f1WMGMdiffmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.18,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.16,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.14,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.12,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.10,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.08,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.18,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.18,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.16,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.16,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.14,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.14,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.12,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.12,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.10,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.10,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.08,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.08,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mode f_1 values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject f_1 mode values','FontSize',18)
    axis([0.5/2 3.5/2 -0.2 0.2])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_f1_mode_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF MD MODE WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==1));
    pCP_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==1));
    pRP_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==1));
    pSC_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==1));
    pSR_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==1));
    pSP_ZOOMint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==1),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==1));
    pCR_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==2));
    pCP_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==2));
    pRP_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==2));
    pSC_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==2));
    pSR_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==2));
    pSP_ZOOMnotint = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==2),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==2));
    pCR_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==3));
    pCP_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==3));
    pRP_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==3));
    pSC_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==1 & MDWMGMdiffmode_grp==3));
    pSR_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==3 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==3));
    pSP_RESOLVE = ranksum(MDWMGMdiffmode_boxplot(subject_grp==4 & MDWMGMdiffmode_grp==3),MDWMGMdiffmode_boxplot(subject_grp==2 & MDWMGMdiffmode_grp==3));
    table_pvals(73,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(78).fig = figure(78);
    set(h(78).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(MDWMGMdiffmode_grp)/2-0.20 unique(MDWMGMdiffmode_grp)/2+0.20]', repmat(MDWMGMdiffmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(MDWMGMdiffmode_grp)/2-0.20 unique(MDWMGMdiffmode_grp)/2+0.20]', repmat(MDWMGMdiffmode_median',1,2)','m-','LineWidth',8)
    scatter(MDWMGMdiffmode_grp(subject_grp==3)/2, MDWMGMdiffmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmode_grp(subject_grp==1)/2, MDWMGMdiffmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmode_grp(subject_grp==4)/2, MDWMGMdiffmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(MDWMGMdiffmode_grp(subject_grp==2)/2, MDWMGMdiffmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.49,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.47,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.45,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.43,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.41,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.39,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.49,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.49,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.47,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.47,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.45,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.45,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.43,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.43,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.41,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.41,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.39,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.39,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mode MD values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject MD mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.50])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_MD_mode_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% VISUALIZATION OF d MODE WM-GM GRADIENT VALUE DISTRIBUTIONS
    pCR_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==1));
    pCP_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==1));
    pRP_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==1));
    pSC_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==1));
    pSR_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==1));
    pSP_ZOOMint = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==1),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==1));
    pCR_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==2));
    pCP_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==2));
    pRP_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==2));
    pSC_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==2));
    pSR_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==2));
    pSP_ZOOMnotint = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==2),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==2));
    pCR_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==3));
    pCP_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==3));
    pRP_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==3));
    pSC_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==1 & dWMGMdiffmode_grp==3));
    pSR_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==3 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==3));
    pSP_RESOLVE = ranksum(dWMGMdiffmode_boxplot(subject_grp==4 & dWMGMdiffmode_grp==3),dWMGMdiffmode_boxplot(subject_grp==2 & dWMGMdiffmode_grp==3));
    table_pvals(74,:) = [pCR_ZOOMint pCP_ZOOMint pRP_ZOOMint pSC_ZOOMint pSR_ZOOMint pSP_ZOOMint ...
        pCR_ZOOMnotint pCP_ZOOMnotint pRP_ZOOMnotint pSC_ZOOMnotint pSR_ZOOMnotint pSP_ZOOMnotint ...
        pCR_RESOLVE pCP_RESOLVE pRP_RESOLVE pSC_RESOLVE pSR_RESOLVE pSP_RESOLVE];
    h(79).fig = figure(79);
    set(h(79).fig, 'Position', [50, 50, 750, 550]);
    plot([unique(dWMGMdiffmode_grp)/2-0.20 unique(dWMGMdiffmode_grp)/2+0.20]', repmat(dWMGMdiffmode_mean',1,2)','c-','LineWidth',8)
    hold on
    plot([unique(dWMGMdiffmode_grp)/2-0.20 unique(dWMGMdiffmode_grp)/2+0.20]', repmat(dWMGMdiffmode_median',1,2)','m-','LineWidth',8)
    scatter(dWMGMdiffmode_grp(subject_grp==3)/2, dWMGMdiffmode_boxplot(subject_grp==3),850, 'k.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmode_grp(subject_grp==1)/2, dWMGMdiffmode_boxplot(subject_grp==1),850, 'r.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmode_grp(subject_grp==4)/2, dWMGMdiffmode_boxplot(subject_grp==4),850, '.', 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerEdgeAlpha',0.7);
    scatter(dWMGMdiffmode_grp(subject_grp==2)/2, dWMGMdiffmode_boxplot(subject_grp==2),850, '.','MarkerEdgeColor',[0 0.749 1], 'jitter','on', 'jitterAmount', 0.16,'MarkerEdgeAlpha',0.7);
    if pCR_ZOOMint < WillThr
        plot((1-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((1-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMint < WillThr
        plot((1-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMint < WillThr
        plot((1-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMint < WillThr
        plot((1-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMint < WillThr
        plot((1-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMint < WillThr
        plot((1-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((1-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((2-0.25)/2,0.83,['CR p = ' num2str(round(pCR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.80,['CM p = ' num2str(round(pCP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pRP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.77,['RM p = ' num2str(round(pRP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSC_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.74,['CS p = ' num2str(round(pSC_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSR_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.71,['RS p = ' num2str(round(pSR_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pSP_ZOOMnotint < WillThr
        plot((2-0.31)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((2-0.25)/2,0.68,['MS p = ' num2str(round(pSP_ZOOMnotint*100000)/100000,4)],'FontSize',14);
    end
    if pCR_RESOLVE < WillThr
        plot((3-0.37)/2,0.83,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15);
        text((3-0.31)/2,0.83,['CR p = ' num2str(round(pCR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pCP_RESOLVE < WillThr
        plot((3-0.37)/2,0.80,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.80,['CM p = ' num2str(round(pCP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pRP_RESOLVE < WillThr
        plot((3-0.37)/2,0.77,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.77,['RM p = ' num2str(round(pRP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSC_RESOLVE < WillThr
        plot((3-0.37)/2,0.74,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.74,['CS p = ' num2str(round(pSC_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSR_RESOLVE < WillThr
        plot((3-0.37)/2,0.71,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.71,['RS p = ' num2str(round(pSR_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    if pSP_RESOLVE < WillThr
        plot((3-0.37)/2,0.68,'k*','LineStyle','none','LineWidth',3,'MarkerSize',15)
        text((3-0.31)/2,0.68,['MS p = ' num2str(round(pSP_RESOLVE*100000)/100000,4)],'FontSize',14);
    end
    hold off
    title({'Distributions of differences between WM and GM';'mode d values from C3-C6 segments over subjects'})
    set(gca,'XTick',(1:1:3)/2,...
         'XTickLabel',{'ZOOMit interp'
                       'ZOOMit non-interp'
                       'RESOLVE'
                       },...
         'TickLength',[0 0],'LineWidth',2,...
         'FontSize',14)
    xlabel('Different protocols','FontSize',18)
    ylabel('Difference of single-subject d mode values','FontSize',18)
    axis([0.5/2 3.5/2 0.0 0.85])
    % axis([0.25 1.75 0.2 0.8])
    grid on
    pause(1)
    print(fullfile(save_path,[num2str(size(subject,1)) '_subjects_d_mode_WM_GM_diff_over_subjects']), '-dpng', '-r300')

    %% SAVE table_pvals VARIABLE INTO dmri_comparison_pvaltable.mat FILE
    close all
    save(fullfile(save_path,'dmri_comparison_pvaltable.mat'))
end
