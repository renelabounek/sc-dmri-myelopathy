function T2TRA_shift_correction_interp( file_name, bindir )
 %--------------------------------------------------------------------
 % Function for slice-by-slice correction T2TRA images from interleaved
 % sequence corrupted by motion artifacts.
 % Steps:
 % 1) Spliting input image into two images - odd and even slices and interpolate empty slices
 % 2) Register images (fixed image - image with original even slices; moving image -image with original odd slices) by isct_antsSliceRegularizedRegistration function (called by slice_reg.sh)
 %    Receive image: odd_to_even
 % 3) Perform combination of registred odd_to_even image and image with fixed even slices
 % Proposed by JV and RL, 2018
 %--------------------------------------------------------------------


 %file_name = '/md2/valosek/2007B_T2TRA/T2TRA_thr_bias_corr'     % uncoment for manual debuging

 % ------------------------------------------------------------------------------------
 % 1) Spliting input image into two images
    if exist([file_name '.nii']) == 0
        disp 'Gunzipping . . .'
        gunzip([file_name '.nii.gz']);
    end

    t2tra = load_untouch_nii([file_name '.nii']);

    [X,Y,Z]=meshgrid(1:size(t2tra.img,1),1:size(t2tra.img,2),1:2:size(t2tra.img,3));
    [Xn,Yn,Zn] = meshgrid(1:size(t2tra.img,1),1:size(t2tra.img,2),2:2:size(t2tra.img,3));

    % Odd fixes, even interp
    data_fixed = t2tra.img(:,:,1:2:size(t2tra.img,3));
    data_interp = interp3(X,Y,Z,data_fixed,Xn,Yn,Zn,'linear');

    data_interp_even = zeros(size(t2tra.img));
    data_interp_even(:,:,1:2:size(t2tra.img,3)) = data_fixed;
    data_interp_even(:,:,2:2:size(t2tra.img,3)) = data_interp;
    data_interp_even(:,:,size(t2tra.img,3)) = t2tra.img(:,:,42);

    % Even fixes, odd interp
    data_fixed = t2tra.img(:,:,2:2:size(t2tra.img,3));
    data_interp = interp3(Xn,Yn,Zn,data_fixed,X,Y,Z,'linear');

    data_interp_odd = zeros(size(t2tra.img));
    data_interp_odd(:,:,1:2:size(t2tra.img,3)) = data_interp;
    data_interp_odd(:,:,2:2:size(t2tra.img,3)) = data_fixed;
    data_interp_odd(:,:,1) = t2tra.img(:,:,1);

    %Save to nii
    t2tra.img = data_interp_even;
    sufix = '_even_interp';
    out=([file_name sufix '.nii']);
    save_untouch_nii(t2tra,out);

    %Save to nii
    t2tra.img = data_interp_odd;
    sufix = '_odd_interp';
    out=([file_name sufix '.nii']);
    save_untouch_nii(t2tra,out);

    %bindir='/md2/NA-CSD/bin';  % uncoment for manual debuging
% ------------------------------------------------------------------------------------
% 2)  Registration by isct_antsSliceRegularizedRegistration
    pathToscript = fullfile(bindir,'slice_reg.sh');
    fixed = ([file_name '_odd_interp.nii']);
    moving = ([file_name '_even_interp.nii']);

    system([pathToscript ' ' fixed ' ' moving]);    % call bash script
% ------------------------------------------------------------------------------------
% 3) Final creation of image by combination registred odd and original even slices

    t2tra = load_untouch_nii([file_name '.nii']);               % necessary load again because was overwritte

    %file_name2 = '/md2/valosek/2007B_T2TRA/odd_to_even'         % uncoment for manual debuging
    file_name2 = fullfile(pwd,'odd_to_even');

    if exist([file_name2 '.nii']) == 0
        disp 'Gunzipping . . .'
        gunzip([file_name2 '.nii.gz']);
    end

    reg = load_untouch_nii([file_name2 '.nii']);

    % Create 3D matrixes containting zeros
    final=zeros(size(t2tra.img));

    % Fill odd slices
    for f=1:2:(size(t2tra.img,3))
        final(:,:,f)=reg.img(:,:,f);
    end

    % Fill even slices
    for f=2:2:(size(t2tra.img,3))
        final(:,:,f)=t2tra.img(:,:,f);
    end

    %Save to nii
    t2tra.img = final;
    sufix = '_reg';
    out=([file_name sufix '.nii']);
    save_untouch_nii(t2tra,out);

    gzip([file_name sufix '.nii']);
    delete([file_name sufix '.nii']);

    delete([file_name '.nii']);

end
