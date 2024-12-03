import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib import pyplot as plt


if __name__ == '__main__':
    ep_obsid = 'ep06800000039wxtCMOS45l23v1'
    obsid = '06800000039'  # 观测号
    cmos = '45'            # 模块号
    s_num = 1              # 源编号
    cts_thres = 10         # 提取图像的阈值，对应的单帧下源区域的光子数目

    cr = 'your data path'
    hdu = fits.open(cr+ep_obsid+'/ep'+obsid+'wxt'+cmos+'.cat')
    hdudata = hdu[1].data
    x0 = hdudata['x'][s_num-1]
    y0 = hdudata['y'][s_num-1]

    x0 = x0 * 16
    y0 = y0 * 16
    hdu = fits.open(cr+ep_obsid+'/ep'+obsid+'wxt'+cmos+'po_cl.evt')
    hdudata = hdu[1].data
    data_cl = np.transpose(np.vstack((hdudata['time'], hdudata['rawx'], hdudata['rawy'], hdudata['x'], hdudata['y'], hdudata['pi'])))
    data_cl = pd.DataFrame(data_cl, columns=['time', 'rawx', 'rawy', 'x', 'y', 'pi'])

    hdu = fits.open(cr+ep_obsid+'/ep'+obsid+'wxt'+cmos+'po_uf.evt')
    hdudata = hdu[1].data
    data_uf = np.transpose(np.vstack((hdudata['time'], hdudata['x'], hdudata['y'], hdudata['pi'], hdudata['cmosfram'])))
    data_uf = pd.DataFrame(data_uf, columns=['time', 'x', 'y', 'pi', 'cmosfram'])

    data = pd.merge(data_cl, data_uf, how='left', on=['time', 'x', 'y', 'pi'])

    fig1 = plt.figure()
    plt.hist2d(data['x'], data['y'], bins=[200, 200])     # 探测图像，可不看，仅供参考

    frame = data['cmosfram'].values
    x = data['x'].values
    y = data['y'].values
    fig2 = plt.figure()
    frame_index = np.arange(np.min(frame), np.max(frame)+1, 1)
    data_hist = np.histogram(frame, bins=frame_index)[0]
    data_hist_all = data_hist
    plt.plot(frame_index[:-1]-frame_index[0], data_hist, color='blue')
#    plt.xlim(0,100000)
    plt.ylim(0,20)
    index = (np.abs(x-x0)<50) & (np.abs(y-y0)<50)
    rawx0 = np.mean(data['rawx'].values[index])
    rawy0 = np.mean(data['rawy'].values[index])
    data_hist = np.histogram(frame[index], bins=frame_index)[0]
    plt.plot(frame_index[:-1]-frame_index[0], data_hist, color='red')
    plt.xlabel('Frame')
    plt.ylabel('Counts')
    # 全帧以及源的帧计数光变曲线，判断源的光变是否有单帧粒子数的增加引起

    # 提取单帧图像，阈值这里设为20，具体参考上一步源光变曲线的计数
    index_all = np.where(data_hist > cts_thres)[0]
    if len(index_all)>0:
        for i in index_all:
            index = (frame == frame_index[i])

            fig3 = plt.figure()
            rawx = data['rawx'][index]
            rawy = data['rawy'][index]
            plt.hist2d(rawx, rawy, bins=[50, 50], range=[[rawx0-200, rawx0+200], [rawy0-200, rawy0+200]])  # 源图像
            fig4 = plt.figure()
            plt.hist2d(rawx, rawy, bins=[100, 100], range=[[0, 4096], [0, 4096]])  # 全帧图像

    plt.show()
