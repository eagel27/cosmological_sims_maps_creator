import numpy as np
import sys

from simulations.IllustrisTNG.read_data import (
    IllustrisTNGGalaxyDataReader,
    IllustrisTNGSimulationReader,
)
from simulations.IllustrisTNG import IllustrisTNGGalaxy
from constants import SIM_NAME, MAPS_RESULTS_PATH, BASE_PATH, CATALOG_DS_PATH
from base.galaxy_base import GalaxyMapsBase


class EagleIllustrisTNGGalaxyMaps(GalaxyMapsBase):
    """
    Custom run of EAGLE simulations analyzed with IllustrisTNG code.
    The data can be loaded with the IllustrisTNG readers.
    """

    SIM = "EAGLE_TNGlike"
    SIM_NAME = SIM_NAME
    BASE_PATH = BASE_PATH
    MAPS_RESULTS_PATH = MAPS_RESULTS_PATH
    CATALOG_DS_PATH = CATALOG_DS_PATH

    def __init__(self, align=False, los="Z", suffix="", id_start=None, id_stop=None):
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
        sim_reader = IllustrisTNGSimulationReader(self.BASE_PATH)
        galaxies_catalog_data = sim_reader.galaxies(snapshot)
        header_params = sim_reader.header_params(snapshot)
        hubble_param = header_params["HubbleParam"]
        subhalo_masses = (
            galaxies_catalog_data["SubhaloMassType"][:, 4] * 1e10 / hubble_param
        )
        subhalo_ids = np.where(subhalo_masses > 1e10)[0]
        subhalo_pos_array = galaxies_catalog_data["SubhaloPos"][subhalo_ids]

        subhalo_data = list(zip(subhalo_ids, subhalo_pos_array))
        return subhalo_data

    def create_galaxy_object(
        self, subhalo_data, snapshot, output_file=None, should_load_datasets=True
    ):
        """
        Create a galaxy object for TNG
        :param subhalo_data: The subhalo data of this galaxy
        :param snapshot: The snapshot that this galaxy exists
        :param output_file: Output file of the maps
        :param should_load_datasets: If this param is True, load all snapshot data
        :return: Return a Galaxy object
        """
        subhalo_id, subhalo_pos = subhalo_data
        galaxy_reader = IllustrisTNGGalaxyDataReader(
            self.BASE_PATH, snapshot, subhalo_id, subhalo_pos
        )
        g_object = IllustrisTNGGalaxy(
            subhalo_id,
            subhalo_pos,
            galaxy_reader,
            align=self.align,
            los=self.los,
            output_file=output_file,
            should_load_datasets=should_load_datasets,
        )
        return g_object
