import numpy as np
import sys

from read_data import IllustrisTNGGalaxyDataReader, IllustrisTNGSimulationReader
from constants import BASE_PATH_FIX, SIM_NAME, MAPS_RESULTS_PATH, CATALOG_DS_PATH
from base.galaxy_base import GalaxyBase, GalaxyMapsBase


class Galaxy(GalaxyBase):
    def __init__(self, subhalo_id, centre, galaxy_snapshot_reader, align=False, los='Z',
                 output_file=None, should_load_datasets=True):
        """
        Initialize all variables for this galaxy
        :param subhalo_id: The subhalo unique ID
        :param centre: The coordinates of the centre of mass
        :param galaxy_snapshot_reader: The IllustrisTNG snapshot reader instance
        """
        self.galaxy_reader = galaxy_snapshot_reader
        self.subhaloID = subhalo_id

        self.hubble_param = self.galaxy_reader.hubble_param
        self.scaling_factor = self.galaxy_reader.scaling_factor
        self.centre = centre * self.scaling_factor / self.hubble_param

        self.len_stars = None
        self.mass_stars = None

        self.align = align
        self.los = los

        self.output_file = output_file
        self.hmrg = self.galaxy_reader.read_catalog_half_mass_radius(4)

        super(Galaxy, self).__init__(should_load_datasets)

    def get_id(self):
        return self.subhaloID

    def get_cosm_params(self):
        cosm_params = {
            'scaling_factor': self.galaxy_reader.scaling_factor,
            'hubble_param': self.galaxy_reader.hubble_param,
            'omega0': self.galaxy_reader.omega0,
            'omegaLambda': self.galaxy_reader.omegaLambda,
            'boxsize': self.galaxy_reader.boxsize
        }
        return cosm_params

    def calc_systemic_velocity(self):
        potential = self.galaxy_potential(4)
        velocities = self.star_velocities

        system_stars_vel = np.mean(velocities[potential < np.percentile(potential, 5)], axis=0)
        system_stars_vel = np.tile(system_stars_vel, [velocities.shape[0], 1])

        return system_stars_vel

    def galaxy_potential(self, itype):
        return self.galaxy_reader.read_galaxy_potential(itype)


class GalaxyMaps(GalaxyMapsBase):

    SIM = 'TNG'
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
        :return: Return a list of tuples (galaxy_id, galaxy_position)
        """
        sim_reader = IllustrisTNGSimulationReader(BASE_PATH_FIX[snapshot].format(snapshot))
        galaxies_catalog_data = sim_reader.galaxies(snapshot)
        header_params = sim_reader.header_params(snapshot)
        hubble_param = header_params['HubbleParam']
        subhalo_masses = galaxies_catalog_data['SubhaloMassType'][:, 4] * 1e10 / hubble_param
        subhalo_ids = np.where(subhalo_masses > 1e10)[0]
        subhalo_pos_array = galaxies_catalog_data['SubhaloPos'][subhalo_ids]

        subhalo_data = list(zip(subhalo_ids, subhalo_pos_array))
        return subhalo_data

    def create_galaxy_object(self, subhalo_data, snapshot, output_file=None, should_load_datasets=True):
        """
        Create a galaxy object for TNG
        :param subhalo_data: The subhalo data of this galaxy
        :param snapshot: The snapshot that this galaxy exists
        :param output_file: Output file of the maps
        :param should_load_datasets: If this param is True, load all snapshot data
        :return: Return a Galaxy object
        """
        subhalo_id, subhalo_pos = subhalo_data
        galaxy_reader = IllustrisTNGGalaxyDataReader(BASE_PATH_FIX[snapshot].format(snapshot),
                                                     snapshot, subhalo_id,
                                                     subhalo_pos)
        g_object = Galaxy(subhalo_id, subhalo_pos, galaxy_reader,
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
