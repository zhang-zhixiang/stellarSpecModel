import os

home_dir = os.path.expanduser('~')
grid_data_dir = '{}/.stellarSpecModel/grid_data/'.format(home_dir)

grid_names = {
    # grid_name: (file_name, url, md5)
    'MARCS': ('MARCS_grid.hdf5', 'https://www.jianguoyun.com/p/DZmcNoUQ2ZfcCBjW-5cFIAA', 'ee32ad18eb0ab427bcebeffb0b4916ce'),
    'BTCond': 'BTCond_grid.hdf5', 'https://www.jianguoyun.com/p/DZmcNoUQ2ZfcCBjW-5cFIAA', 'ee32ad18eb0ab427bcebeffb0b4916ce'),
}