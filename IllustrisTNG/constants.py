import os.path

DEBUG = False

SNAPSHOT = 99
SIM_NAME = 'TNG100-1'

if DEBUG:
    BASE = '/home/eirini/Documents/PhD'
    BASE_PATH = '/home/eirini/Documents/PhD/Data/IllustrisTNG/{}/output'.format(SIM_NAME)
    BASE_PATH_FIX = {
        99: '/home/eirini/Documents/PhD/Data/IllustrisTNG/{}/output'.format(SIM_NAME),
        91: '/home/eirini/Documents/PhD/Data/IllustrisTNG/{}/output'.format(SIM_NAME)
    }

else:
    BASE = '/net/diva/scratch1/eirinia/projects'
    BASE_PATH = '/net/diva/scratch1/anegri/simulations/TNG/%s/0{}' % (SIM_NAME, )
    BASE_PATH_FIX = {
        99: '/storage/scratch/anegri/IllustrisTNG/%s/0{}' % (SIM_NAME, ),
        91: '/storage/scratch/eirinia/TNG/%s/0{}' % (SIM_NAME, )
    }


CATALOG_DS_PATH = '/u/aeirini/Data/TNG100/postprocessing/balanced_dataset/balanced_alignments.hdf5'
MAPS_RESULTS_PATH = os.path.join(BASE, 'Results/Datasets/Maps/IllustrisTNG/{}/'.format(SIM_NAME))

if SIM_NAME.startswith('Illustris'):
    little_h = 0.704
else:
    little_h = 0.6774


