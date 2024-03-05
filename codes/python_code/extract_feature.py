import os
from multiprocessing import Pool
import scipy.io as scio
import numpy as np
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy import signal
import h5py


def get_time_data(data_ori, is_interpolate=True):
    # print(data_ori)
    if is_interpolate:
        x = np.linspace(1, 10, 10)
        func = interpolate.UnivariateSpline(x, data_ori, s=0)  # 强制通过所有点
        x = np.linspace(1, 10, 100)
        data_ori = func(x)
    maxI = max(data_ori)
    minI = min(data_ori)
    tMAXI = list(data_ori).index(maxI)
    spline = UnivariateSpline(x, data_ori - maxI / 2, s=0)
    roots = spline.roots()  # find the roots
    W = abs(roots[1] - roots[0]) if len(roots) >= 2 else 1e6  # 如果没有半高宽则说明曲线单调递增或递减，即半高宽无限大
    L = maxI / W
    numP = len(signal.argrelextrema(data_ori, np.greater)[0])
    numT = len(signal.argrelextrema(data_ori, np.less)[0])
    std = np.std(data_ori, ddof=1)
    var = np.var(data_ori)
    S = 0
    K = 0
    for i in range(len(data_ori)):
        S += (data_ori[i] - std) ** 3
        K += (data_ori[i] - std) ** 4
    S = S / (len(data_ori) * (var ** 3))
    K = K / (len(data_ori) * (var ** 4)) - 3
    # if w == 0 or var == 0:
    #     plt.figure()
    #     x = list(range(len(data_ori)))
    #     plt.plot(x, data_ori)
    #     plt.show()
    return [maxI, minI, tMAXI, W, L, numP, numT, S, K]


def get_all_data(data, mask, mean_data, step):
    Data = []
    (row, col, sli, frame) = data.shape
    for i in range(1, row - 1, step):
        for j in range(1, col - 1, step):
            for k in range(1, sli - 1, step):
                # location and label
                data_list = [i, j, k, mask[i, j, k]]
                # space features
                # 6 领域
                nei_6 = [mean_data[i, j, k], mean_data[i - 1, j, k], mean_data[i + 1, j, k], mean_data[i, j - 1, k],
                         mean_data[i, j + 1, k], mean_data[i, j, k - 1], mean_data[i, j, k + 1]]
                mean_value = np.mean(nei_6)
                std_value = np.std(nei_6)
                min_value = np.min(nei_6)
                max_value = np.max(nei_6)
                data_list += [mean_value, std_value, abs(mean_data[i, j, k] - min_value),
                              abs(mean_data[i, j, k] - max_value)]
                # 26 邻域
                nei_26 = mean_data[i - 1:i + 2, j - 1:j + 2, k - 1:k + 2]
                mean_value = np.mean(nei_26)
                std_value = np.std(nei_26)
                min_value = np.min(nei_26)
                max_value = np.max(nei_26)
                data_list += [mean_value, std_value, abs(mean_data[i, j, k] - min_value),
                              abs(mean_data[i, j, k] - max_value)]

                # time features
                data_ori = data[i, j, k, :]
                if np.var(data_ori) != 0:
                    Data.append(data_list + get_time_data(data_ori))
    return Data


def main_fun(file_path, mask_path, feature_path, mean_img):
    print(file_path, ' start')
    # scio read image
    data = scio.loadmat(file_path)
    data = data['iSNAP_dMRA_tmp200']
    # h5py read image
    # data = h5py.File(file_path)
    # data = np.array(data['iSNAP_dMRA_tmp200'])
    # data = data.transpose((3, 2, 1, 0))  # h5py read image need to update shape
    print('image shape:', data.shape)
    mask = scio.loadmat(mask_path)
    mask = mask['final_seg']
    print('mask shape:', data.shape)
    mean_data = np.mean(data, axis=3)
    print('mip image shape:', mean_data.shape)
    scio.savemat(mean_img, {'mean_data': mean_data})
    Data = get_all_data(data, mask, mean_data, 1)
    np.savetxt(feature_path, Data, fmt='%f', delimiter=',')
    print(file_path, ' end')


if __name__ == '__main__':
    dir_1 = 'D:\Reconstructed_iSNAP_images_MMD'
    dir_2 = 'D:\Reconstructed_iSNAP_images_Exp20220423'

    case1 = ['MMD20210513_LYH', 'MMD20210621_GMP', 'MMD20210707_PostSurg_LBY', 'MMD20210707_PostSurg_YCY']
    case2 = os.listdir(dir_2)

    img_filename = 'iSNAP_dMRA_tmp200.mat'
    # 请注意，该mask并不是真实的分割或者标注mask
    # 该mask仅用于对voxel进行分层抽样
    # 可以使用阈值等简单算法获得该mask
    mask_filename = 'new_dMRA200_vISMRM_final_seg.mat'

    # file_dir = [os.path.join(dir_1, i) for i in case1] + [os.path.join(dir_2, i, 'RAW', 'vessels') for i in case2]
    file_dir = [os.path.join(dir_1, i) for i in case1]
    # file_dir = ['D:\一分区\Reconstructed_iSNAP_images_MMD\MMD20210707_PostSurg_LBY']
    p = Pool(5)
    for case_dir in file_dir:
        file_path = os.path.join(case_dir, img_filename)
        mask_path = os.path.join(case_dir, mask_filename)
        feature_path = os.path.join(case_dir, 'last_features.txt')
        mean_img = os.path.join(case_dir, 'last_mean_img.mat')
        # main_fun(file_path, mask_path, feature_path, mean_img)
        p.apply_async(main_fun, args=(file_path, mask_path, feature_path, mean_img))
    p.close()
    p.join()
