%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE: vessel segmentation
% chanyo@ohsu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc
close all;


% put data into /DATA/
smpl_name = dir('./DATA/ROI_*');
marker_name = {'D240','SMA','LYVE1','CD34','AQP1','DAPI','Mask'};
    
output_path = './Result/';
mkdir(output_path);

for i=1:length(smpl_name)
   
    fname = dir(sprintf('./DATA/%s/*.tif', smpl_name(i).name));
    
    % read compositive figure
    sub_ = dir(sprintf('./DATA/%s/*Other*/*composite*', smpl_name(i).name));
    Icomp = imread( sprintf('%s/%s', sub_.folder, sub_.name));
    
    for k=1:3 % adjust intensity
        Icomp(:,:,k) = imadjust(Icomp(:,:,k));
    end
    
    % image file name
    for j=1:length(fname)
        fname_{j} = sprintf('%s',fname(j).name);
    end
    
    fig_order = [];
    for j=1:length(marker_name)
        id = []; id = find ( contains( fname_, marker_name{j}, 'IgnoreCase',true) > 0);
        
        I(:,:,j) = imread( sprintf('./DATA/%s/%s', smpl_name(i).name, fname(id).name));
        
        if strcmp(marker_name{j}, 'DAPI')
            mask_DAPI = imbinarize(imgaussfilt(I(:,:,j),2));
        end
    end
   
    
    
    I_z =  uint8(zeros(size(I(:,:,1))));
    % compare LYVE1 and CD34
    if length(find (imbinarize(I(:,:,3))> 0)) > length(find (imbinarize(I(:,:,4))> 0))
        marker_idx = [1 2 4];
    else
        marker_idx = [1 2 3];
    end
    
    % projected image
    for ii=1:length(marker_idx)
        j = marker_idx(ii);
        I_ = gray2ind((I(:,:,j)));
        I_ = imadjust(I_);
        
        I_z = I_z +  imgaussfilt(I_,3);
    end
    
    % segmentation
    bw = imadjust(I_z);
    bw = imbinarize(bw);
    bw = imfill(bw,'holes');
    mask = bwareafilt(bw, [300 inf]);
    mask = imdilate(mask,strel('disk',1));

    
   
    
    
    %% remove DAPI seg
    mask_f = mask;
    s = regionprops('table',mask, mask_DAPI,'MeanIntensity','PixelIdxList');
    cutoff_DAPI = 0.75; % cutoff for DAPI overlap
    
    id = []; id = find(s.MeanIntensity >= cutoff_DAPI);
    
    for j=1:length(id)
        mask_f(s.PixelIdxList{ id(j)}) = 0;
    end
   
    
    L = bwlabeln(mask_f);
   
    %imagesc(mask); hold on, visboundaries(mask_f);
    
    
    %% feature extraction
    T = [];
    T = regionprops('table',L, 'Area','ConvexArea','FilledArea','MajoraxisLength', 'MinoraxisLength',...
                    'EquivDiameter','Perimeter', 'Extent','Eccentricity','Centroid');
    s = []; s = regionprops('table',L,L,'MeanIntensity');             
    Vessel_ID = s.MeanIntensity;
    T.VesselID = Vessel_ID;
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,1)),'MeanIntensity');             
    s = []; s = regionprops('table',L,(I(:,:,1)),'MeanIntensity');             
    Int_D240 = s.MeanIntensity;
    T.Int_D240 = Int_D240;
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,2)),'MeanIntensity');             
    s = []; s = regionprops('table',L,(I(:,:,2)),'MeanIntensity');             
    Int_aSMA = s.MeanIntensity;
    T.Int_aSMA = Int_aSMA;
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,3)),'MeanIntensity');             
    s = []; s = regionprops('table',L,(I(:,:,3)),'MeanIntensity');             
    Int_LYVE1 = s.MeanIntensity;
    T.Int_LYVE1 = Int_LYVE1;
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,4)),'MeanIntensity');    
    s = []; s = regionprops('table',L,(I(:,:,4)),'MeanIntensity');             
    Int_CD34 = s.MeanIntensity;
    T.Int_CD34 = Int_CD34;
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,5)),'MeanIntensity');             
    s = []; s = regionprops('table',L,(I(:,:,5)),'MeanIntensity');             
    Int_AQP1 = s.MeanIntensity;
    T.Int_AQP1 = Int_AQP1;
    
    
    %s = []; s = regionprops('table',L,gray2ind(I(:,:,6)),'MeanIntensity');             
    s = []; s = regionprops('table',L,(I(:,:,6)),'MeanIntensity');             
    Int_DAPI = s.MeanIntensity;
    T.Int_DAPI = Int_DAPI;
    
    writetable(T, sprintf('%sFeat_%s.csv', output_path, smpl_name(i).name));
    
    %% result image
    
    figure('pos',[10 10 1000 600]);
    for j=1:length(fname)
        ax(j)=subplot(2,4,j);
        
        imagesc(imadjust(gray2ind(I(:,:,j)))); hold on, 
        visboundaries(mask,'Color','g'); % segmented mask

        visboundaries(I(:,:,end)); % provided mask
        title(sprintf('%s',marker_name{j}));
        colormap gray

    end

    linkaxes(ax);
    
%     figure; 
%     imagesc(Icomp); hold on; visboundaries(mask,'Color','w');
%     visboundaries(I(:,:,6),'Color','y');
%     
    
    overlay_mask = imoverlay(Icomp, bwperim(mask_f==1), [1 1 1]);
    overlay_mask = imoverlay(overlay_mask, bwperim(I(:,:,end)==1), [1 0 0]);
    imwrite(overlay_mask, sprintf('%soverlaid_mask_%s.tif', output_path, smpl_name(i).name), 'tif');

    
    imwrite(mask_f, sprintf('%svessel_mask_%s.tif',output_path, smpl_name(i).name),'tif');
    imwrite(uint8(L), sprintf('%svessel_mask_label_%s.tif',output_path, smpl_name(i).name),'tif');

end

