function fill_holes( image )

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
