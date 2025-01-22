import os

DB_USER = ""
DB_PASSWORD = ""
SIM_NAME = 'EAGLE_TNGlike'
DEBUG = True

SNAPSHOT = 28

RESULTS_PATH = '/home/eirini/Documents/PhD/Results/'
if DEBUG:
    BASE_PATH = '/home/eirini/Documents/PhD'
    DATA_PATH = os.path.join(BASE_PATH, 'Data/EAGLE/{}'.format(SIM_NAME))
    MAPS_RESULTS_PATH = os.path.join(BASE_PATH, 'Results/Datasets/Maps/EAGLE_TNGlike/')
    SIM_PATH = os.path.join(BASE_PATH, 'Data/EAGLE/')
    NFILES = 16
else:
    BASE_PATH = '/net/diva/scratch1/eirinia/'
    DATA_PATH = os.path.join(BASE_PATH, 'projects/Data/EAGLE/{}'.format(SIM_NAME))
    MAPS_RESULTS_PATH = os.path.join(BASE_PATH, 'projects/Results/Datasets/Maps/EAGLE/{}/'.format(SIM_NAME))
    SIM_PATH = os.path.join(BASE_PATH, 'projects/Data/EAGLE/')
    NFILES = 256


STARS_CUTOUT_Z0 = os.path.join(DATA_PATH, 'Extracted/stars_cutout_snap_028.h5py')
INSITU_STARS = os.path.join(DATA_PATH, 'DB_Data/insitu_all.h5py')

SNAPSHOT_PATH_SUFFIX = {
    28: 'snapshot_028_z000p000/snap_028_z000p000',
    27: 'snapshot_027_z000p101/snap_027_z000p101',
    26: 'snapshot_026_z000p183/snap_026_z000p183',
    25: 'snapshot_025_z000p271/snap_025_z000p271'
}

SNAPSHOT_PATH = os.path.join(DATA_PATH, '%s.{}.hdf5')
SNAPSHOT_PATH_Z0 = SNAPSHOT_PATH % SNAPSHOT_PATH_SUFFIX[SNAPSHOT]

CATALOG_DS_PATH = '/u/aeirini/Data/EAGLE/postprocessing/balanced_dataset/balanced_alignments.hdf5'
