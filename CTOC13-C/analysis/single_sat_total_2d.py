# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

'''
Analysis of single satellite data
- Data:
    Shape: (10000000, 7)
    Columns: [sma, ecc, inc, raan, argp, ta, n_visible]
    Special: ecc = argp = ta = 0
- Analysis:
    1. sma: certain height
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

if __name__ == '__main__':
    # Save data
    data_name = 'single_sat_1108_2_large'
    # data_bin_path = os.path.join(os.path.dirname(__file__), '../data', data_name + '.bin')
    # data = read_data(data_bin_path)
    # np.save(os.path.join(os.path.dirname(__file__), '../data', data_name + '.npy'), data)

    # Load data
    data_np_path = os.path.join(os.path.dirname(__file__), '../data', data_name + '.npy')
    data = np.load(data_np_path)
    print(data.shape)

    # Analysis
    step = 125
    sma_idx = 40
    sma = data[sma_idx:data.shape[0]:step, 0]
    inc = np.rad2deg(data[sma_idx:data.shape[0]:step, 2])
    raan = np.rad2deg(data[sma_idx:data.shape[0]:step, 3])
    n_visible = data[sma_idx:data.shape[0]:step, 6]

    # 1. inc & raan vs n_visible: 2D heatmap
    plt.figure()
    heatmap, xedges, yedges = np.histogram2d(inc, raan, bins=100, weights=n_visible)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.imshow(heatmap.T, extent=extent, origin='lower', cmap='jet', aspect='auto')
    cbar = plt.colorbar()
    cbar.set_label('n_visible')

    # Adjust colorbar ticks to reflect original n_visible values
    tick_locs = np.linspace(heatmap.min(), heatmap.max(), num=10)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels([f'{int(val)}' for val in np.linspace(n_visible.min(), n_visible.max(), num=10)])

    plt.xlabel('inc (deg)')
    plt.ylabel('raan (deg)')
    plt.title('inc & raan vs n_visible, sma = {}'.format(sma[0]))

    plt.show()