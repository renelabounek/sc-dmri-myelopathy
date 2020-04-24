%% SIGNATURE
% Implemented by Rene Labounek
% fMRI Laboratory, Department of Neurology, Palacky University and University Hospital Olomouc, Czech Republic
% Division of Clinical Behavioral Neuroscience, Department of Pediatrics, University of Minnesota, Minneapolis, USA
% contact emails: rlaboune@umn.edu, rene.labounek@gmail.com
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