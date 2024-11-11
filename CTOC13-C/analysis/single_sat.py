# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm

'''
Analysis of single satellite data
- Data:
    Shape: (10000000, 7)
    Columns: [sma, ecc, inc, raan, argp, ta, n_visible]
    Special: ecc = argp = ta = 0
- Analysis:
    1. sma: 6878.137 < sma < 7378.137
    2. inc: 0.0 < inc < PI
    3. raan: 0.0 < raan < 2 * PI
- Plot:
    1. sma & inc & raan vs n_visible
    2. inc vs n_visible
    3. raan vs n_visible
'''


# Read data from binary file
def read_data(db_filepath):
    with open(db_filepath, 'rb') as file:
        binary_data = file.read()
        text_data = binary_data.decode('utf-8').replace('\n', ' ').replace('\t', ' ')
        text_data = text_data.strip()
        data = np.array(text_data.split(), dtype=np.float64)
    data = data.reshape(-1, 7)
    return data

def custom_colormap():
    cmap = (ListedColormap(['blue', 'royalblue', 'lightgreen', 'yellow', 'orange'])
        .with_extremes(over='red'))
    
    return cmap


if __name__ == '__main__':
    # Save data
    data_name = 'single_sat_ship'
    # data_bin_path = os.path.join(os.path.dirname(__file__), '../data', data_name + '.bin')
    # data = read_data(data_bin_path)
    # np.save(os.path.join(os.path.dirname(__file__), '../data', data_name + '.npy'), data)

    # Load data
    data_np_path = os.path.join(os.path.dirname(__file__), '../data', data_name + '.npy')
    data = np.load(data_np_path)
    print(data.shape)

    # Analysis
    step = 100
    sma = data[0:data.shape[0]:step, 0]
    inc = np.rad2deg(data[0:data.shape[0]:step, 2])
    raan = np.rad2deg(data[0:data.shape[0]:step, 3])
    n_visible = data[0:data.shape[0]:step, 6]

    # 1. sma & inc & raan vs n_visible: 3D heatmap
    # Point color: n_visible
    # Color bar shown on the right
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    # 创建自定义颜色映射
    custom_cmap = custom_colormap()
    
    # 使用自定义的颜色映射和归一化
    bounds = np.linspace(0, 8, 5)
    norm = BoundaryNorm(bounds, custom_cmap.N, extend='max')
    
    scatter = ax.scatter(sma, inc, raan, c=n_visible, cmap=custom_cmap, norm=norm)
    
    ax.set_xlabel('sma(km)')
    ax.set_ylabel('inc(deg)')
    ax.set_zlabel('raan(deg)')

    fig.colorbar(scatter)
    plt.title('sma & inc & raan vs n_ship_visible (max 8)')
   
    plt.show()


