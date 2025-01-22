import numpy as np
import os
import sys
import copy
import pandas as pd

sys.path.append(os.path.abspath('../..'))

from db_queries import QUERY_SUBHALO_DATA
from read_data import EagleReader, EagleSnapshotDataReader, EagleGalaxyDataReader
from constants import DB_USER, DB_PASSWORD, SIM_NAME, MAPS_RESULTS_PATH, DATA_PATH, CATALOG_DS_PATH
from base.galaxy_base import GalaxyBase, GalaxyMapsBase


class Galaxy(GalaxyBase):

    def __init__(self, gid, gn, sgn, centre, centre_vel, hmrg, galaxy_snapshot_reader,
                 los='Z', align=False, output_file=None, should_load_datasets=True):
        """
        Initialize all variables for this galaxy
        :param gid: The GalaxyID
        :param gn: The GroupID
        :param sgn: The SubGroupID
        :param centre: The coordinates of the centre of mass
        :param centre_vel: The velocities (3D) of centre of potential
        :param hmrg: The half-mass radius of this galaxy
        :param galaxy_snapshot_reader: The EAGLE snapshot reader instance
        """
        self.gid = gid
        self.gn = gn
        self.sgn = sgn
        self.centre = centre
        self.centre_vel = centre_vel
        self.hmrg = hmrg

        self.len_stars = None
        self.mass_stars = None

        self.align = align
        self.los = los
        self.galaxy_reader = galaxy_snapshot_reader
        self.output_file = output_file

        self.SIMNAME = 'EAGLE'
        super(Galaxy, self).__init__(should_load_datasets)

    def get_id(self):
        return self.gid

    def get_cosm_params(self):
        cosm_params = {
            'scaling_factor': self.galaxy_reader.scaling_factor,
            'hubble_param': self.galaxy_reader.hubble_param,
            'omega0': self.galaxy_reader.omega0,
            'omegaLambda': self.galaxy_reader.omegaLambda,
            'boxsize': self.galaxy_reader.boxsize
        }
        return cosm_params

    def calc_systemic_velocity(self, percentage=0.1):
        coords = self.star_coordinates
        velocities = self.star_velocities

        central_mask = np.linalg.norm(coords, axis=1) <= percentage * self.hmrg
        vel_central = velocities[central_mask]
        vel_central = np.mean(vel_central, axis=0)

        return vel_central


class GalaxyMaps(GalaxyMapsBase):

    SIM = 'EAGLE'
    SIM_NAME = SIM_NAME
    MAPS_RESULTS_PATH = MAPS_RESULTS_PATH
    CATALOG_DS_PATH = CATALOG_DS_PATH

    def __init__(self, align=False, los='Z', suffix='', id_start=None, id_stop=None):
        """
        Initialize params required to create 2D maps
        :param align: If align is True, create the edge on view of the galaxy
        :param los: If the LOS is specified, use it for the maps creation
        """
        self.align = align
        self.los = los
        self.suffix = suffix
        self.id_start = id_start
        self.id_stop = id_stop

    def load_catalog_data(self, snapshot):
        """
        Load the catalog data of the galaxies
        :return: Return the query data for that galaxy (in list form)
        """
        sim_reader = EagleReader(DB_USER, DB_PASSWORD, snapshot)
        galaxies_catalog_data_path = os.path.join(DATA_PATH,
                                                  'DB_Data/galaxy_catalog_data_{}.npy'.format(snapshot))
        if not os.path.exists(galaxies_catalog_data_path):
            galaxies_catalog_data = sim_reader.query(SIM_NAME, QUERY_SUBHALO_DATA, extra=(snapshot, ))
            np.save(galaxies_catalog_data_path, galaxies_catalog_data)
        else:
            galaxies_catalog_data = np.load(galaxies_catalog_data_path)

        # sim_reader.snapshot_reader.save_cutout()
        del sim_reader
        return galaxies_catalog_data

    def create_galaxy_object(self, subhalo_data, snapshot,
                             output_file=None, should_load_datasets=True):
        """
        Create a galaxy object for EAGLE
        :param subhalo_data: The subhalo data of this galaxy
        :param output_file: Output file of the maps
        :param should_load_datasets: If this param is True, load all snapshot data
        :return: Return a Galaxy object
        """
        snapshot_reader = EagleSnapshotDataReader(snapshot)
        snapshot_reader.read_header()

        centre = subhalo_data[7:10]
        # convert from Mpc to kpc and physical
        centre = np.array(centre) * snapshot_reader.a * 1e3
        centre_vel = subhalo_data[10:13]
        hmrg = subhalo_data[14]

        galaxy_reader = None
        if should_load_datasets:
            galaxy_reader = EagleGalaxyDataReader(*subhalo_data[0:3], centre,
                                                  centre_vel, snapshot_reader)
        g_object = Galaxy(*subhalo_data[0:3], centre, centre_vel, hmrg, galaxy_reader,
                          align=self.align, los=self.los, output_file=output_file,
                          should_load_datasets=should_load_datasets)

        return g_object


if __name__ == '__main__':
    mode = int(sys.argv[1])
    align = int(sys.argv[2])

    id_start, id_stop = None, None
    suffix = ''
    if len(sys.argv) > 3:
        id_start = int(sys.argv[3])
        id_stop = int(sys.argv[4])
        suffix = sys.argv[5]

    if align == 1:
        galaxy_maps = GalaxyMaps(align=True, los='X', suffix=suffix, id_start=id_start, id_stop=id_stop)
    else:
        galaxy_maps = GalaxyMaps(suffix=suffix, id_start=id_start, id_stop=id_stop)

    if mode == 1:
        galaxy_maps.create_maps()
    elif mode == 2:
        galaxy_maps.save_2d_maps()
