from sklearn.covariance import EllipticEnvelope
from sklearn.model_selection import GridSearchCV, train_test_split
import warnings
import pickle
from sklearn.ensemble import IsolationForest
from sklearn.metrics import make_scorer, recall_score
import os
import scipy.io as scio
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
import numpy as np
import time
import random

random.seed(123)
warnings.filterwarnings("ignore")


def predict_data(case_name, model, location, data, seg_data, model_name, feature):
    print(case_name, 'start')
    predict_mask = np.zeros(seg_data.shape)
    Data_pre = model.predict(np.nan_to_num(data.astype(np.float32)))
    for i in range(location.shape[0]):
        if Data_pre[i] == -1:
            predict_mask[int(location[i, 0]), int(location[i, 1]), int(location[i, 2])] = 1
    mask_save_path = '../save_mask/' + model_name + '/' + feature + '/'
    if not os.path.exists(mask_save_path):
        os.makedirs(mask_save_path)
    scio.savemat(os.path.join(mask_save_path, case_name + '.mat'), {'predict_mask': predict_mask})
    print(case_name, ' end')


def std_deal(data):
    for i in range(data.shape[1]):
        data[:, i] = (data[:, i] - data[:, i].min()) / (data[:, i].max() - data[:, i].min())
    return data


def extract_data(feature_path, mask_path, feature, is_train, sample_config=None):
    Data = np.loadtxt(feature_path, delimiter=',')
    seg_data = scio.loadmat(mask_path)['final_seg']

    Data[np.isnan(Data)] = 0
    Data[np.isinf(Data)] = 0

    label = Data[:, 3]
    if is_train:
        if sample_config['down_sample']:
            if not sample_config['balance']:
                _, data, _, label = train_test_split(Data, label, test_size=sample_config['sample_rate'], random_state=123, stratify=label)
            else:
                vessel_Data = Data[label == 1, :]
                vessel_len = vessel_Data.shape[0]
                background_Data = Data[label == 0, :]
                np.random.seed(123)
                row_rand_array = np.arange(background_Data.shape[0])
                np.random.shuffle(row_rand_array)
                background_Data_select = background_Data[row_rand_array[:vessel_len]]
                data = np.vstack((vessel_Data, background_Data_select))
                label = data[:, 3]
    else:
        data = Data
    location = data[:, :3]
    if feature == 'time':
        data = data[:, 12:]
    elif feature == 'space':
        data = data[:, 4:12]
    elif feature == 'time_space':
        data = data[:, 4:]
    return location, data, label, seg_data


if __name__ == '__main__':
    model_name = 'IsolationForest'
    dir_1 = 'D:\Reconstructed_iSNAP_images_MMD'
    dir_2 = 'D:\Reconstructed_iSNAP_images_Exp20220423'

    case1 = ['MMD20210513_LYH', 'MMD20210621_GMP', 'MMD20210707_PostSurg_LBY', 'MMD20210707_PostSurg_YCY']
    case2 = os.listdir(dir_2)
    case1.sort()
    case2.sort()

    file_dir_1 = [os.path.join(dir_1, i) for i in case1]
    file_dir_2 = [os.path.join(dir_2, i) for i in case2]

    feature_filename = 'last_features.txt'
    mask_filename = 'new_dMRA200_vISMRM_final_seg.mat'

    for file_dirs in [file_dir_1, file_dir_2]:
        train_feature_path = os.path.join(file_dirs[0], feature_filename)
        train_mask_path = os.path.join(file_dirs[0], mask_filename)
        sample_config = {'down_sample': True, 'balance': False, 'sample_rate': 0.01}
        feature = 'time_space'
        _, data, label, _ = extract_data(train_feature_path, train_mask_path, feature,
                                         is_train=True, sample_config=sample_config)
        contamination = label.sum() / len(label)
        label[label == 0] = -1
        recall_fraud = make_scorer(recall_score, pos_label=-1)

        gs_params = {
            # 'max_samples': [300, 500, 1000],
            'contamination': [contamination],
            # 'max_features': [1, 3, 7],
            # 'n_estimators': [1000],
            'random_state': [1]
        }
        clf = IsolationForest()

        start_time = time.time()
        grid = GridSearchCV(clf, gs_params, verbose=1, cv=5, scoring=recall_fraud)
        grid.fit(data, label)

        print(grid.best_params_)
        print(grid.best_score_)
        subject = 'MMD' if 'MMD' in file_dirs[0].split('\\')[2] else 'NOR'
        save_dir = '../save_model/' + feature + '/' + subject + '/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        save_model_file = save_dir + model_name + '_' + str(grid.best_score_) + '.pkl'
        best_model = grid.best_estimator_
        with open(save_model_file, 'wb') as f:
            pickle.dump(best_model, f)

        with open(save_model_file, 'rb') as f:
            model = pickle.load(f)
        model.fit(data)
        print(model_name, 'end train, time cost:', time.time() - start_time)

        for file_dir in file_dirs:
            train_feature_path = os.path.join(file_dir, feature_filename)
            train_mask_path = os.path.join(file_dir, mask_filename)
            location, data, _, seg_data = extract_data(train_feature_path, train_mask_path, feature, is_train=False)
            predict_data(file_dir.split('\\')[2], model, location, data, seg_data, model_name, feature)
