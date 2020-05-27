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
clear all;close all;clc;
%% Data read
demographic_file='/home/user/subject_table.xlsx';
[num, txt, raw] = xlsread(demographic_file);
age=[raw{2:end,4}]';
grp=[raw{2:end,2}]';
gender= {raw{2:end,3}}';

%% Subgroup stats
F1=0;
F2=0;
F3=0;
F4=0;

for ind = 1:size(gender,1)
        if grp(ind,1)==1 && strcmp(gender{ind,1},'F')
                F1 = F1 + 1;
        elseif grp(ind,1)==2 && strcmp(gender{ind,1},'F')
                F2 = F2 + 1;
        elseif grp(ind,1)==3 && strcmp(gender{ind,1},'F')
                F3 = F3 + 1;
        elseif grp(ind,1)==4 && strcmp(gender{ind,1},'F')
                F4 = F4 + 1;
        end
end
numADCCCF = F2+F4;
numADCCC = sum(grp==2 | grp == 4);
numHC = sum(grp==1);
numHCF = F1;
numADCCCMF = F2;
numADCCCM = sum(grp==2);
numADCCCSF = F4;
numADCCCS = sum(grp == 4);
numREPRODUCIBILITYF = F3/2;
numREPRODUCIBILITY = sum(grp == 3)/227.4;
clear F1 F2 F3 F4

ageHC=age(grp==1);
ageHCstat=[mean(ageHC) std(ageHC)];
ageADCCC=age(grp==2 | grp==4);
ageADCCCstat=[mean(ageADCCC) std(ageADCCC)];
ageADCCCM=age(grp==2);
ageADCCCMstat=[mean(ageADCCCM) std(ageADCCCM)];
ageADCCCS=age(grp==4);
ageADCCCSstat=[mean(ageADCCCS) std(ageADCCCS)];
ageREPRODUCIBILITY=age(grp==3);
ageREPRODUCIBILITYstat=[mean(ageREPRODUCIBILITY) std(ageREPRODUCIBILITY)];

p(1,1) = ranksum(ageHC,ageHC);
[~, p(1,2)] = ttest2(ageHC,ageHC);
p(2,1) = ranksum(ageHC,ageADCCC);
[~, p(2,2)] = ttest2(ageHC,ageADCCC);
p(3,1) = ranksum(ageHC,ageADCCCM);
[~, p(3,2)] = ttest2(ageHC,ageADCCCM);
p(4,1) = ranksum(ageHC,ageADCCCS);
[~, p(4,2)] = ttest2(ageHC,ageADCCCS);
p(5,1) = ranksum(ageHC,ageREPRODUCIBILITY);
[~, p(5,2)] = ttest2(ageHC,ageREPRODUCIBILITY);
