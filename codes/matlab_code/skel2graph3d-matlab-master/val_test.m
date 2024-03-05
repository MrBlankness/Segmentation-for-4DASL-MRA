% clc;
% clear;
% path = 'res_DT_position_001.csv';
% [poi1,~]=importdata(path);
% path = 'res_RF_position_001.csv';
% [poi2,~]=importdata(path);
% path = 'res_float_position_001.csv';
% [poi3,~]=importdata(path);
% 
maskPath = 'D:\idmWorkspace\Data Set 1 Digital subtraction angiography images and aneurysm mask and geometry data for tasks 1, 2, and 3_2\CADA-Training_MaskImages-NIFTI\A001_masks.nii.gz';
nii=load_untouch_nii(maskPath);
maskData=nii.img;
figure('Name','真实肿瘤点'),volshow(maskData);

maskPath = 'D:\idmWorkspace\Data Set 1 Digital subtraction angiography images and aneurysm mask and geometry data for tasks 1, 2, and 3_2\CADA-Training_MaskImages-NIFTI\A001_labeledMasks.nii.gz';
nii=load_untouch_nii(maskPath);
maskData=nii.img;
figure('Name','真实肿瘤点2'),volshow(maskData);

% 
% 
% poiData = maskData * 0;
% for i = 1:size(poi1,1)
%     poiData(poi1(i,1),poi1(i,2),poi1(i,3))=1;
% end
% figure('Name','DT'),volshow(poiData);
% res = poiData.*maskData;
% sum(sum(sum(res)))
% 
% poiData = maskData * 0;
% for i = 1:size(poi2,1)
%     poiData(poi2(i,1),poi2(i,2),poi2(i,3))=1;
% end
% figure('Name','RF'),volshow(poiData);
% res = poiData.*maskData;
% sum(sum(sum(res)))
% 
% poiData = maskData * 0;
% for i = 1:size(poi3,1)
%     poiData(poi3(i,1),poi3(i,2),poi3(i,3))=1;
% end
% figure('Name','float'),volshow(poiData);
% res = poiData.*maskData;
% sum(sum(sum(res)))