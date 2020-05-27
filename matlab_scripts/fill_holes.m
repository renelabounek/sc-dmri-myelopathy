function fill_holes( image )
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
    T2TRA_bin = load_untouch_nii([image '.nii']);
    %imshow4(T2TRA_bin.img);

    for i=1:size(T2TRA_bin.img,3)
        T2TRA_bin_fill(:,:,i) = imfill(T2TRA_bin.img(:,:,i),'holes');
    end

    for i=1:size(T2TRA_bin.img,3)
        T2TRA_bin_fill_remove(:,:,i) = bwareaopen(T2TRA_bin_fill(:,:,i),30000);
    end

    % imshow4(T2TRA_bin_fill_remove);

    T2TRA_bin.img = T2TRA_bin_fill_remove;
    out=[image '_fill.nii'];
    save_untouch_nii(T2TRA_bin,out);

end
