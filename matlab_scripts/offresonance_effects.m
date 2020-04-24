%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
clc;
clear all;
close all;
%% Result import
save_path='/home/user/results';
workspace_file = fullfile(save_path,'dmri_comparison_pvaltable.mat');
load(workspace_file);
%% Mean and STD value estimations
offres_c3c6(1,1) = mean(fieldcoef_median(:,1));
offres_c3c6(2,1)  = std(fieldcoef_median(:,1));

temp = fieldcoef_median(:,2);
temp(temp==0) = [];
offres_c3c6(1,2) = mean(temp);
offres_c3c6(2,2)  = std(temp);

offres_c3c6(1,3) = mean(fieldcoef_median(:,3));
offres_c3c6(2,3)  = std(fieldcoef_median(:,3));

offres_c3c6 = round(offres_c3c6)

offres_c3(1,1) = mean(fieldcoef3_median(:,1));
offres_c3(2,1)  = std(fieldcoef3_median(:,1));

temp = fieldcoef3_median(:,2);
temp(temp==0) = [];
offres_c3(1,2) = mean(temp);
offres_c3(2,2)  = std(temp);

offres_c3(1,3) = mean(fieldcoef3_median(:,3));
offres_c3(2,3)  = std(fieldcoef3_median(:,3));

offres_c3 = round(offres_c3)

offres_c5c6(1,1) = mean(fieldcoef56_median(:,1));
offres_c5c6(2,1)  = std(fieldcoef56_median(:,1));

temp = fieldcoef56_median(:,2);
temp(temp==0) = [];
offres_c5c6(1,2) = mean(temp);
offres_c5c6(2,2)  = std(temp);

offres_c5c6(1,3) = mean(fieldcoef56_median(:,3));
offres_c5c6(2,3)  = std(fieldcoef56_median(:,3));

offres_c5c6 = round(offres_c5c6)