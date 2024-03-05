clc;
clear;
close all;

addpath('./skel2graph3d-matlab-master');
addpath('./skeleton3d-matlab-master');

data_dir = 'predict_mask/';
case_names = dir([data_dir, '*.mat']);
for case_i = 1:length(case_names)
    case_name = case_names(case_i).name;
    path = [data_dir, case_name];
    load(path);
    predict_mask = getMaxVol(predict_mask);
    skel = Skeleton3D(imbinarize(double(predict_mask)));
    length = cal_skel_length(skel);
    fprintf('%.2f, ', length * 0.8); % voxel size 0.8 mm isotropic
end
fprintf('\n');
