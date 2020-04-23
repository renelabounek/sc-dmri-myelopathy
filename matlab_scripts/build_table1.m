%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com

save_path='/home/user/results';
workspace_file = fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);

subnum=60;

group_table(1,1) = median(sqrt(FA_wm_var([subject{:,2}]==1,1)));
group_table(1,2) = median(sqrt(FA_wm_var([subject{:,2}]==2,1)));
group_table(1,4) = table_pvals(9,2);
group_table(1,5) = median(sqrt(FA_wm_var([subject{:,2}]==4,1)));
group_table(1,7) = table_pvals(9,4);
group_table(1,8) = median(sqrt(FA_wm_var([subject{:,2}]==1,3)));
group_table(1,9) = median(sqrt(FA_wm_var([subject{:,2}]==2,3)));
group_table(1,11) = table_pvals(9,14);
group_table(1,12) = median(sqrt(FA_wm_var([subject{:,2}]==4,3)));
group_table(1,14) = table_pvals(9,16);

group_table(2,1) = median(FA_wm_kurtosis([subject{:,2}]==1,1));
group_table(2,2) = median(FA_wm_kurtosis([subject{:,2}]==2,1));
group_table(2,4) = table_pvals(63,2);
group_table(2,5) = median(FA_wm_kurtosis([subject{:,2}]==4,1));
group_table(2,7) = table_pvals(63,4);
group_table(2,8) = median(FA_wm_kurtosis([subject{:,2}]==1,3));
group_table(2,9) = median(FA_wm_kurtosis([subject{:,2}]==2,3));
group_table(2,11) = table_pvals(63,14);
group_table(2,12) = median(FA_wm_kurtosis([subject{:,2}]==4,3));
group_table(2,14) = table_pvals(63,16);

group_table(3,1) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==1,1)));
group_table(3,2) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==2,1)));
group_table(3,4) = table_pvals(28,2);
group_table(3,5) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==4,1)));
group_table(3,7) = table_pvals(28,4);
group_table(3,8) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==1,3)));
group_table(3,9) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==2,3)));
group_table(3,11) = table_pvals(28,14);
group_table(3,12) = median(sqrt(f1_gauss_wm_var([subject{:,2}]==4,3)));
group_table(3,14) = table_pvals(28,16);

group_table(4,1) = median(f1_gauss_wm_skewness([subject{:,2}]==1,1));
group_table(4,2) = median(f1_gauss_wm_skewness([subject{:,2}]==2,1));
group_table(4,4) = table_pvals(57,2);
group_table(4,5) = median(f1_gauss_wm_skewness([subject{:,2}]==4,1));
group_table(4,7) = table_pvals(57,4);
group_table(4,8) = median(f1_gauss_wm_skewness([subject{:,2}]==1,3));
group_table(4,9) = median(f1_gauss_wm_skewness([subject{:,2}]==2,3));
group_table(4,11) = table_pvals(57,14);
group_table(4,12) = median(f1_gauss_wm_skewness([subject{:,2}]==4,3));
group_table(4,14) = table_pvals(57,16);

group_table(5,1) = median(f1_gauss_wm_kurtosis([subject{:,2}]==1,1));
group_table(5,2) = median(f1_gauss_wm_kurtosis([subject{:,2}]==2,1));
group_table(5,4) = table_pvals(65,2);
group_table(5,5) = median(f1_gauss_wm_kurtosis([subject{:,2}]==4,1));
group_table(5,7) = table_pvals(65,4);
group_table(5,8) = median(f1_gauss_wm_kurtosis([subject{:,2}]==1,3));
group_table(5,9) = median(f1_gauss_wm_kurtosis([subject{:,2}]==2,3));
group_table(5,11) = table_pvals(65,14);
group_table(5,12) = median(f1_gauss_wm_kurtosis([subject{:,2}]==4,3));
group_table(5,14) = table_pvals(65,16);

ptblind=20;
group_table(6,1) = median(MD_gauss_wm_median([subject{:,2}]==1,1));
group_table(6,2) = median(MD_gauss_wm_median([subject{:,2}]==2,1));
group_table(6,4) = table_pvals(ptblind,2);
group_table(6,5) = median(MD_gauss_wm_median([subject{:,2}]==4,1));
group_table(6,7) = table_pvals(ptblind,4);
group_table(6,8) = median(MD_gauss_wm_median([subject{:,2}]==1,3));
group_table(6,9) = median(MD_gauss_wm_median([subject{:,2}]==2,3));
group_table(6,11) = table_pvals(ptblind,14);
group_table(6,12) = median(MD_gauss_wm_median([subject{:,2}]==4,3));
group_table(6,14) = table_pvals(ptblind,16);

ptblind=45;
group_table(7,1) = median(MD_gauss_wm_mean([subject{:,2}]==1,1));
group_table(7,2) = median(MD_gauss_wm_mean([subject{:,2}]==2,1));
group_table(7,4) = table_pvals(ptblind,2);
group_table(7,5) = median(MD_gauss_wm_mean([subject{:,2}]==4,1));
group_table(7,7) = table_pvals(ptblind,4);
group_table(7,8) = median(MD_gauss_wm_mean([subject{:,2}]==1,3));
group_table(7,9) = median(MD_gauss_wm_mean([subject{:,2}]==2,3));
group_table(7,11) = table_pvals(ptblind,14);
group_table(7,12) = median(MD_gauss_wm_mean([subject{:,2}]==4,3));
group_table(7,14) = table_pvals(ptblind,16);

ptblind=23;
group_table(8,1) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==1,1)));
group_table(8,2) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==2,1)));
group_table(8,4) = table_pvals(ptblind,2);
group_table(8,5) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==4,1)));
group_table(8,7) = table_pvals(ptblind,4);
group_table(8,8) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==1,3)));
group_table(8,9) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==2,3)));
group_table(8,11) = table_pvals(ptblind,14);
group_table(8,12) = median(sqrt(MD_gauss_wm_var([subject{:,2}]==4,3)));
group_table(8,14) = table_pvals(ptblind,16);

ptblind=33;
group_table(9,1) = median(sqrt(d_gauss_wm_var([subject{:,2}]==1,1)));
group_table(9,2) = median(sqrt(d_gauss_wm_var([subject{:,2}]==2,1)));
group_table(9,4) = table_pvals(ptblind,2);
group_table(9,5) = median(sqrt(d_gauss_wm_var([subject{:,2}]==4,1)));
group_table(9,7) = table_pvals(ptblind,4);
group_table(9,8) = median(sqrt(d_gauss_wm_var([subject{:,2}]==1,3)));
group_table(9,9) = median(sqrt(d_gauss_wm_var([subject{:,2}]==2,3)));
group_table(9,11) = table_pvals(ptblind,14);
group_table(9,12) = median(sqrt(d_gauss_wm_var([subject{:,2}]==4,3)));
group_table(9,14) = table_pvals(ptblind,16);

ptblind=61;
group_table(10,1) = median(d_gauss_wm_skewness([subject{:,2}]==1,1));
group_table(10,2) = median(d_gauss_wm_skewness([subject{:,2}]==2,1));
group_table(10,4) = table_pvals(ptblind,2);
group_table(10,5) = median(d_gauss_wm_skewness([subject{:,2}]==4,1));
group_table(10,7) = table_pvals(ptblind,4);
group_table(10,8) = median(d_gauss_wm_skewness([subject{:,2}]==1,3));
group_table(10,9) = median(d_gauss_wm_skewness([subject{:,2}]==2,3));
group_table(10,11) = table_pvals(ptblind,14);
group_table(10,12) = median(d_gauss_wm_skewness([subject{:,2}]==4,3));
group_table(10,14) = table_pvals(ptblind,16);

ptblind=21;
group_table(11,1) = median(MD_gauss_gm_median([subject{:,2}]==1,1));
group_table(11,2) = median(MD_gauss_gm_median([subject{:,2}]==2,1));
group_table(11,4) = table_pvals(ptblind,2);
group_table(11,5) = median(MD_gauss_gm_median([subject{:,2}]==4,1));
group_table(11,7) = table_pvals(ptblind,4);
group_table(11,8) = median(MD_gauss_gm_median([subject{:,2}]==1,3));
group_table(11,9) = median(MD_gauss_gm_median([subject{:,2}]==2,3));
group_table(11,11) = table_pvals(ptblind,14);
group_table(11,12) = median(MD_gauss_gm_median([subject{:,2}]==4,3));
group_table(11,14) = table_pvals(ptblind,16);

ptblind=46;
group_table(12,1) = median(MD_gauss_gm_mean([subject{:,2}]==1,1));
group_table(12,2) = median(MD_gauss_gm_mean([subject{:,2}]==2,1));
group_table(12,4) = table_pvals(ptblind,2);
group_table(12,5) = median(MD_gauss_gm_mean([subject{:,2}]==4,1));
group_table(12,7) = table_pvals(ptblind,4);
group_table(12,8) = median(MD_gauss_gm_mean([subject{:,2}]==1,3));
group_table(12,9) = median(MD_gauss_gm_mean([subject{:,2}]==2,3));
group_table(12,11) = table_pvals(ptblind,14);
group_table(12,12) = median(MD_gauss_gm_mean([subject{:,2}]==4,3));
group_table(12,14) = table_pvals(ptblind,16);

ptblind=60;
group_table(13,1) = median(MD_gauss_gm_skewness([subject{:,2}]==1,1));
group_table(13,2) = median(MD_gauss_gm_skewness([subject{:,2}]==2,1));
group_table(13,4) = table_pvals(ptblind,2);
group_table(13,5) = median(MD_gauss_gm_skewness([subject{:,2}]==4,1));
group_table(13,7) = table_pvals(ptblind,4);
group_table(13,8) = median(MD_gauss_gm_skewness([subject{:,2}]==1,3));
group_table(13,9) = median(MD_gauss_gm_skewness([subject{:,2}]==2,3));
group_table(13,11) = table_pvals(ptblind,14);
group_table(13,12) = median(MD_gauss_gm_skewness([subject{:,2}]==4,3));
group_table(13,14) = table_pvals(ptblind,16);

ptblind=31;
group_table(14,1) = median(d_gauss_gm_median([subject{:,2}]==1,1));
group_table(14,2) = median(d_gauss_gm_median([subject{:,2}]==2,1));
group_table(14,4) = table_pvals(ptblind,2);
group_table(14,5) = median(d_gauss_gm_median([subject{:,2}]==4,1));
group_table(14,7) = table_pvals(ptblind,4);
group_table(14,8) = median(d_gauss_gm_median([subject{:,2}]==1,3));
group_table(14,9) = median(d_gauss_gm_median([subject{:,2}]==2,3));
group_table(14,11) = table_pvals(ptblind,14);
group_table(14,12) = median(d_gauss_gm_median([subject{:,2}]==4,3));
group_table(14,14) = table_pvals(ptblind,16);

ptblind=51;
group_table(15,1) = median(d_gauss_gm_mean([subject{:,2}]==1,1));
group_table(15,2) = median(d_gauss_gm_mean([subject{:,2}]==2,1));
group_table(15,4) = table_pvals(ptblind,2);
group_table(15,5) = median(d_gauss_gm_mean([subject{:,2}]==4,1));
group_table(15,7) = table_pvals(ptblind,4);
group_table(15,8) = median(d_gauss_gm_mean([subject{:,2}]==1,3));
group_table(15,9) = median(d_gauss_gm_mean([subject{:,2}]==2,3));
group_table(15,11) = table_pvals(ptblind,14);
group_table(15,12) = median(d_gauss_gm_mean([subject{:,2}]==4,3));
group_table(15,14) = table_pvals(ptblind,16);

ptblind=22;
group_table(16,1) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==1,1));
group_table(16,2) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==2,1));
group_table(16,4) = table_pvals(ptblind,2);
group_table(16,5) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==4,1));
group_table(16,7) = table_pvals(ptblind,4);
group_table(16,8) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==1,3));
group_table(16,9) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==2,3));
group_table(16,11) = table_pvals(ptblind,14);
group_table(16,12) = median(MD_gauss_diff_wm_gm_median([subject{:,2}]==4,3));
group_table(16,14) = table_pvals(ptblind,16);

ptblind=49;
group_table(17,1) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==1,1));
group_table(17,2) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==2,1));
group_table(17,4) = table_pvals(ptblind,2);
group_table(17,5) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==4,1));
group_table(17,7) = table_pvals(ptblind,4);
group_table(17,8) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==1,3));
group_table(17,9) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==2,3));
group_table(17,11) = table_pvals(ptblind,14);
group_table(17,12) = median(MD_gauss_diff_wm_gm_mean([subject{:,2}]==4,3));
group_table(17,14) = table_pvals(ptblind,16);

ptblind=32;
group_table(18,1) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==1,1));
group_table(18,2) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==2,1));
group_table(18,4) = table_pvals(ptblind,2);
group_table(18,5) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==4,1));
group_table(18,7) = table_pvals(ptblind,4);
group_table(18,8) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==1,3));
group_table(18,9) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==2,3));
group_table(18,11) = table_pvals(ptblind,14);
group_table(18,12) = median(d_gauss_diff_wm_gm_median([subject{:,2}]==4,3));
group_table(18,14) = table_pvals(ptblind,16);

ptblind=52;
group_table(19,1) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==1,1));
group_table(19,2) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==2,1));
group_table(19,4) = table_pvals(ptblind,2);
group_table(19,5) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==4,1));
group_table(19,7) = table_pvals(ptblind,4);
group_table(19,8) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==1,3));
group_table(19,9) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==2,3));
group_table(19,11) = table_pvals(ptblind,14);
group_table(19,12) = median(d_gauss_diff_wm_gm_mean([subject{:,2}]==4,3));
group_table(19,14) = table_pvals(ptblind,16);

ptblind=75;
pom=heu_FAWM(1:subnum);
pom(:,3)=heu_FAWM(end-subnum+1:end);
group_table(20,1) = median(pom([subject{:,2}]==1,1));
group_table(20,2) = median(pom([subject{:,2}]==2,1));
group_table(20,4) = table_pvals(ptblind,2);
group_table(20,5) = median(pom([subject{:,2}]==4,1));
group_table(20,7) = table_pvals(ptblind,4);
group_table(20,8) = median(pom([subject{:,2}]==1,3));
group_table(20,9) = median(pom([subject{:,2}]==2,3));
group_table(20,11) = table_pvals(ptblind,14);
group_table(20,12) = median(pom([subject{:,2}]==4,3));
group_table(20,14) = table_pvals(ptblind,16);

ptblind=76;
pom=heu_f1WM(1:subnum);
pom(:,3)=heu_f1WM(end-subnum+1:end);
group_table(21,1) = median(pom([subject{:,2}]==1,1));
group_table(21,2) = median(pom([subject{:,2}]==2,1));
group_table(21,4) = table_pvals(ptblind,2);
group_table(21,5) = median(pom([subject{:,2}]==4,1));
group_table(21,7) = table_pvals(ptblind,4);
group_table(21,8) = median(pom([subject{:,2}]==1,3));
group_table(21,9) = median(pom([subject{:,2}]==2,3));
group_table(21,11) = table_pvals(ptblind,14);
group_table(21,12) = median(pom([subject{:,2}]==4,3));
group_table(21,14) = table_pvals(ptblind,16);

ptblind=77;
pom=heu_MDWM(1:subnum);
pom(:,3)=heu_MDWM(end-subnum+1:end);
group_table(22,1) = median(pom([subject{:,2}]==1,1));
group_table(22,2) = median(pom([subject{:,2}]==2,1));
group_table(22,4) = table_pvals(ptblind,2);
group_table(22,5) = median(pom([subject{:,2}]==4,1));
group_table(22,7) = table_pvals(ptblind,4);
group_table(22,8) = median(pom([subject{:,2}]==1,3));
group_table(22,9) = median(pom([subject{:,2}]==2,3));
group_table(22,11) = table_pvals(ptblind,14);
group_table(22,12) = median(pom([subject{:,2}]==4,3));
group_table(22,14) = table_pvals(ptblind,16);

ptblind=78;
pom=heu_dWM(1:subnum);
pom(:,3)=heu_dWM(end-subnum+1:end);
group_table(23,1) = median(pom([subject{:,2}]==1,1));
group_table(23,2) = median(pom([subject{:,2}]==2,1));
group_table(23,4) = table_pvals(ptblind,2);
group_table(23,5) = median(pom([subject{:,2}]==4,1));
group_table(23,7) = table_pvals(ptblind,4);
group_table(23,8) = median(pom([subject{:,2}]==1,3));
group_table(23,9) = median(pom([subject{:,2}]==2,3));
group_table(23,11) = table_pvals(ptblind,14);
group_table(23,12) = median(pom([subject{:,2}]==4,3));
group_table(23,14) = table_pvals(ptblind,16);

ptblind=81;
pom=heu_MDGM(1:subnum);
pom(:,3)=heu_MDGM(end-subnum+1:end);
group_table(24,1) = median(pom([subject{:,2}]==1,1));
group_table(24,2) = median(pom([subject{:,2}]==2,1));
group_table(24,4) = table_pvals(ptblind,2);
group_table(24,5) = median(pom([subject{:,2}]==4,1));
group_table(24,7) = table_pvals(ptblind,4);
group_table(24,8) = median(pom([subject{:,2}]==1,3));
group_table(24,9) = median(pom([subject{:,2}]==2,3));
group_table(24,11) = table_pvals(ptblind,14);
group_table(24,12) = median(pom([subject{:,2}]==4,3));
group_table(24,14) = table_pvals(ptblind,16);

ptblind=82;
pom=heu_dGM(1:subnum);
pom(:,3)=heu_dGM(end-subnum+1:end);
group_table(25,1) = median(pom([subject{:,2}]==1,1));
group_table(25,2) = median(pom([subject{:,2}]==2,1));
group_table(25,4) = table_pvals(ptblind,2);
group_table(25,5) = median(pom([subject{:,2}]==4,1));
group_table(25,7) = table_pvals(ptblind,4);
group_table(25,8) = median(pom([subject{:,2}]==1,3));
group_table(25,9) = median(pom([subject{:,2}]==2,3));
group_table(25,11) = table_pvals(ptblind,14);
group_table(25,12) = median(pom([subject{:,2}]==4,3));
group_table(25,14) = table_pvals(ptblind,16);

group_table(:,3) = (group_table(:,2)-group_table(:,1) )./group_table(:,1)*100;
group_table(:,6) = (group_table(:,5)-group_table(:,1) )./group_table(:,1)*100;
group_table(:,10) = (group_table(:,9)-group_table(:,8) )./group_table(:,8)*100;
group_table(:,13) = (group_table(:,12)-group_table(:,8) )./group_table(:,8)*100;