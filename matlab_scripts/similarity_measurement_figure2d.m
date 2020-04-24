%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
clc;
clear all;
close all;
%% INPUT PARAMETER SETTINGS
data_folder='/home/user/data';
% data_folder='/md6';
save_path='/home/user/results';
read_data = 1;
plot_results = 1;
WillThr = 0.05/6;

workspace_file = fullfile(save_path,'dmri_similarity_256_bins_wmgm.mat');
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
        'S0018'	1;...
        'S0019'	2;...
        'S0020'	1;...
        'S0021'	1;...
        'S0023'	1;...
        'S0024'	3;...
        'S0026'	1;...
        'S0027'	3;...
        'S0028'	1;...
        'S0029'	2;...
        'S0030'	2;...
        'S0033'	1;...
        'S0034'	1;...
        'S0035'	2;...
        'S0036'	2;...
        'S0037'	2;...
        'S0038'	2;...
        'S0039'	2;...
        'S0040'	2;...
        'S0041'	2;...
        'S0042'	2;...
        'S0043'	1;...
        'S0044'	2;...
        'S0045'	2;...
        'S0046'	1;...
        'S0047'	1;...
        'S0048'	4;...
        'S0049'	4;...
        'S0050'	4;...
        'S0051'	4;...
        'S0052'	4;...
        'S0053'	4;...
        'S0054'	4;...
        'S0056'	4;...
        'S0057'	4;...
        'S0058'	4;...
        'S0059'	4;...
        'S0060'	2;...
        'S0061'	4;...
        'S0062'	4;...
        'S0063'	3;...
        'S0064'	3;...
        'S0065'	3;...
        'S0066'	3;...
        'S0067'	3 ...
        };
    sbj = cell2mat(subject(:,2));

    protocol = {'ZOOMit_interp'; 'ZOOMit_notinterp' ;'RESOLVE'};
    name={'zoomit'; 'zoomit'; 'resolve'};

%% Similarity estimation
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

                wmgm = 2*mask + gm;
                wmgm(wmgm==3) = 1;

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

                    d_file = fullfile(dmri_preproc_folder,'diff_in_t2tra','mean_dsamples.nii');
                    gunzip([d_file '.gz']);
                    d_file_hdr = spm_vol(d_file);
                    d = spm_read_vols(d_file_hdr);
                    delete(d_file)

                    FA_vec = FA(mask==1);
%                         T2TRA_vec = T2TRA(mask==1);
                    T2TRA_vec = wmgm(mask==1);
                    f1_vec = f1(mask==1);
                    MD_vec = MD(mask==1);
                    MD_vec_norm = (MD_vec - min(MD_vec)) / (max(MD_vec)-min(MD_vec));
                    d_vec = d(mask==1);
                    d_vec_norm = (d_vec - min(d_vec)) / (max(d_vec)-min(d_vec));


                    cor = corrcoef(FA_vec,T2TRA_vec);
                    cor_protocols(sub,ind) = cor(1,2);

                    n=size(FA_vec,1);
                    t_value(sub,ind) = cor_protocols(sub,ind) / sqrt( (1-cor_protocols(sub,ind)^2)/(n-2)  );
                    p_value(sub,ind) = 2*(1-tcdf(abs(t_value(sub,ind)),n-1));

                    T2TRA_vec_norm = round(((T2TRA_vec-min(T2TRA_vec)) / ( max(T2TRA_vec) - min(T2TRA_vec)))*255)+1;
                    FA_vec_norm = round(((FA_vec-min(FA_vec)) / ( max(FA_vec) - min(FA_vec)))*255)+1;
                    MD_vec_norm = round(((MD_vec-min(MD_vec)) / ( max(MD_vec) - min(MD_vec)))*255)+1;
                    f1_vec_norm = round(((f1_vec-min(f1_vec)) / ( max(f1_vec) - min(f1_vec)))*255)+1;
                    d_vec_norm = round(((d_vec-min(d_vec)) / ( max(d_vec) - min(d_vec)))*255)+1;

                    mi_FA(sub,ind) = mi(FA_vec_norm,T2TRA_vec_norm,256);
                    mi_MD(sub,ind) = mi(MD_vec_norm,T2TRA_vec_norm,256);
                    mi_d(sub,ind) = mi(d_vec_norm,T2TRA_vec_norm,256);
                    mi_f1(sub,ind) = mi(f1_vec_norm,T2TRA_vec_norm,256);
                end
                disp(['Subject ' num2str(sub) ' session ' num2str(ind) ' done.'])
            end                
    end
    save(workspace_file,'cor_protocols','t_value','p_value','mi_FA','mi_MD','mi_d','mi_f1')
end
%% Similarity visualization
if plot_results == 1
    load(workspace_file)
    mi_d_boxplot = mi_d(:);
    grp = repmat(1:size(mi_d,2),size(mi_d,1),1);
    grp = grp(:);
    grp(mi_d_boxplot==0) = [];
    mi_d_boxplot(mi_d_boxplot==0) = [];
    mi_MD_boxplot = mi_MD(:);
    mi_MD_boxplot(mi_MD_boxplot==0) = [];
    mi_FA_boxplot = mi_FA(:);
    mi_FA_boxplot(mi_FA_boxplot==0) = [];
    mi_f1_boxplot = mi_f1(:);
    mi_f1_boxplot(mi_f1_boxplot==0) = [];

    h(1).fig = figure(1);
    set(h(1).fig, 'Position', [50, 50, 1200, 280]);
    figure(1)
    subplot(1,4,4)
    bh = boxplot(mi_d_boxplot,grp);
    for i=1:size(bh,1) % <- # graphics handles/x
        for j = 1:size(bh,2)
            set(bh(i,j),'linewidth',2);
            pause(.1);
        end
    end
    title('d')
    grid on
    axis([0.5 3.5 0 0.35])
    set(gca,'XTick',1:size(mi_d,2),...
            'XTickLabel',{'HZi' 'HZni' 'DR'},...
            'TickLength',[0 0],...
            'FontSize',14,...
            'LineWidth',2)
    subplot(1,4,3)
    bh = boxplot(mi_MD_boxplot,grp);
    for i=1:size(bh,1) % <- # graphics handles/x
        for j = 1:size(bh,2)
            set(bh(i,j),'linewidth',2);
            pause(.1);
        end
    end
    title('MD')
    grid on
    axis([0.5 3.5 0 0.35])
    set(gca,'XTick',1:size(mi_d,2),...
            'XTickLabel',{'HZi' 'HZni' 'DR'},...
            'TickLength',[0 0],...
            'FontSize',14,...
            'LineWidth',2)
    subplot(1,4,1)
    bh = boxplot(mi_FA_boxplot,grp);
    for i=1:size(bh,1) % <- # graphics handles/x
        for j = 1:size(bh,2)
            set(bh(i,j),'linewidth',2);
            pause(.1);
        end
    end
    title('FA')
%         ylabel('Mutual Information')
    grid on
    axis([0.5 3.5 0 0.15])
    set(gca,'XTick',1:size(mi_d,2),...
            'XTickLabel',{'HZi' 'HZni' 'DR'},...
            'TickLength',[0 0],...
            'FontSize',14,...
            'LineWidth',2)
    subplot(1,4,2)
    bh = boxplot(mi_f1_boxplot,grp);
    for i=1:size(bh,1) % <- # graphics handles/x
        for j = 1:size(bh,2)
            set(bh(i,j),'linewidth',2);
            pause(.1);
        end
    end
    title('f1')
    grid on
    axis([0.5 3.5 0 0.15])
    set(gca,'XTick',1:size(mi_d,2),...
            'XTickLabel',{'HZi' 'HZni' 'DR'},...
            'TickLength',[0 0],...
            'FontSize',14,...
            'LineWidth',2)
end