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
