import illustris_python as ill
import numpy as np


class IllustrisTNGGalaxyDataReader(object):
    """
    Reads data for a specific galaxy from the IllustrisTNG simulation given the base path.
    The galaxy is identified by the subhalo_id and the snapshot.
    """
    def __init__(self, basePath, snapshot, subhalo_id, subhalo_pos):
        self.basePath = basePath
        self.snapshot = snapshot
        self.subhaloID = subhalo_id

        self.centre = subhalo_pos
        self.header = self.load_header(snapshot)

        self.catalog_fields = None
        self.stars_mask = None
        self.scaling_factor = self.read_header_value('Time')
        self.omega0 = self.read_header_value('Omega0')
        self.omegaLambda = self.read_header_value('OmegaLambda')
        self.hubble_param = self.read_header_value('HubbleParam')
        self.boxsize = self.read_header_value('BoxSize') * self.scaling_factor / self.hubble_param

        # convert to physical kpc
        self.centre *= self.scaling_factor / self.hubble_param

    def load_header(self, snapshot):
        header = ill.groupcat.loadHeader(self.basePath, snapshot)
        return header

    def read_header_value(self, value):
        return self.header.get(value, 1)

    def read_galaxy_from_catalog(self):
        self.catalog_fields = ill.groupcat.loadSingle(self.basePath,
                                                      self.snapshot,
                                                      subhaloID=self.subhaloID)

    def read_catalog_galaxy_mass(self, itype):
        if self.catalog_fields is None:
            self.read_galaxy_from_catalog()
        return self.catalog_fields['SubhaloMassType'][itype]

    def read_catalog_half_mass_radius(self, itype):
        if self.catalog_fields is None:
            self.read_galaxy_from_catalog()
        return self.catalog_fields['SubhaloHalfmassRadType'][itype] * self.scaling_factor / self.hubble_param

    def read_catalog_galaxy_nparticles(self, itype):
        if self.catalog_fields is None:
            self.read_galaxy_from_catalog()
        return self.catalog_fields['SubhaloLenType'][itype]

    def read_catalog_galaxy_fields(self, fields):
        if self.catalog_fields is None:
            self.read_galaxy_from_catalog()
        return [self.catalog_fields[key] for key in fields]

    def read_galaxy_dm_mass(self):
        # Create an array of length n_particles each set to dm_mass.
        mass_total_dm = self.read_catalog_galaxy_mass(1)
        n_particles_dm = self.read_catalog_galaxy_nparticles(1)
        dm_mass = mass_total_dm / n_particles_dm
        dens = np.ones(n_particles_dm, dtype='f8') * dm_mass
        return dens

    def read_galaxy_coordinates(self, part_type, wrap_around_centre=True):
        """ Resolve coordinates of type of particles for this galaxy.
        Coordinates are wrapped around the centre to account for periodicity, if required. """
        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot,
                                        self.subhaloID,  part_type, fields=["Coordinates",
                                                                            "GFM_StellarFormationTime"])
        coords = data['Coordinates'] * self.scaling_factor / self.hubble_param

        # Periodic wrap coordinates around centre.
        if wrap_around_centre:
            boxsize = self.boxsize
            coords = np.mod(coords - self.centre + 0.5 * boxsize, boxsize) + self.centre - 0.5 * boxsize

        if part_type == 4:
            self.stars_mask = np.where(data['GFM_StellarFormationTime'] > 0)[0]
            coords = coords[self.stars_mask]

        # Relative to centre
        coords -= self.centre
        return coords

    def read_galaxy_particle_ids(self, part_type):
        """ Resolve particle ids of part_type particles for this galaxy. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        part_type, fields=["ParticleIDs"])
        if part_type == 4:
            #self.stars_mask = np.where(data['GFM_StellarFormationTime'] > 0)[0]
            data = data[self.stars_mask]
        return data

    def read_galaxy_velocities(self, part_type):
        """ Resolve velocities of part_type particles for this galaxy. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        part_type, fields=["Velocities"])

        data *= np.sqrt(self.scaling_factor)

        if part_type == 4:
            data = data[self.stars_mask]

        return data

    def read_galaxy_potential(self, part_type):
        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        part_type, fields=['Potential'])
        data = data[self.stars_mask]
        return data

    def read_galaxy_mass(self, part_type):
        """ Resolve mass of part_type particles for this galaxy. """
        if part_type != 1:
            data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                            part_type, fields=["Masses"]) * 1e10 / self.hubble_param
            if part_type == 4:
                data = data[self.stars_mask]
        else:
            data = self.read_galaxy_dm_mass()
        return data

    def read_galaxy_metallicity(self):
        """ Resolve metallicity of star particles for this galaxy. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        'stars', fields=["GFM_Metallicity"])
        data = data[self.stars_mask]

        return data

    def read_galaxy_stellar_formation_time(self):
        """ Resolve formation time of star particles for this galaxy
        (expansion factor when this star was born). """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        'stars', fields=["GFM_StellarFormationTime"])
        data = data[self.stars_mask]
        return data

    def read_galaxy_gas_density(self):
        """ Resolve density of gas particles for this galaxy. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        'gas', fields=["Density"])
        return data

    def read_galaxy_sfr(self):
        """ Resolve star formation rate of each gas particle for this galaxy. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        'gas', fields=["StarFormationRate"])
        return data

    def read_galaxy_fields(self, part_type, fields):
        """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        extract the requested attribute of all particles of a selected type. """

        data = ill.snapshot.loadSubhalo(self.basePath, self.snapshot, self.subhaloID,
                                        part_type, fields=fields)
        return data


class IllustrisTNGSnapshotDataReader(object):
    """
    Reads data for a specific snapshot for the IllustrisTNG simulation given the base path.
    """
    def __init__(self, path):
        self.path = path

        # Initialize header variables
        self.a, self.h, self.boxsize = None, None, None
        self.omega0, self.omegaLambda, self.omegaBaryon = None, None, None

        self.dataset_dict = {}

    def read_dataset(self, snapshot, part_type, fields):
        """ Read a selected subset of the data """
        data = ill.snapshot.loadSubset(self.path, snapshot, part_type, fields)
        return data

    def read_subhalo(self, snapshot, subhalo_id, part_type, fields):
        """ Read a selected subhalo (specific fields) """
        data = ill.snapshot.loadSubhalo(self.path, snapshot, subhalo_id, part_type, fields)
        return data

    def read_halo(self, snapshot, halo_id, part_type, fields):
        """ Read a selected halo (specific fields) """
        data = ill.snapshot.loadHalo(self.path, snapshot, halo_id, part_type, fields)
        return data


class IllustrisTNGCatalogReader(object):
    """
    Reads catalog data for the IllustrisTNG simulation given the base path.
    """
    def __init__(self, base_path):
        self.basePath = base_path

    def load_header(self, snapshot):
        header = ill.groupcat.loadHeader(self.basePath, snapshot)
        return header

    def load_halos(self, snapshot, fields):
        halos = ill.groupcat.loadHalos(self.basePath, snapshot, fields=fields)
        return halos

    def load_subhalos(self, snapshot, fields=None):
        subhalos = ill.groupcat.loadSubhalos(self.basePath, snapshot, fields=fields)
        return subhalos


class IllustrisTNGSimulationReader(object):
    """
    Reads data for IllustrisTNG simulation given the base path.
    """
    def __init__(self, base_path):
        self.catalog_reader = IllustrisTNGCatalogReader(base_path)
        self.snapshot_reader = IllustrisTNGSnapshotDataReader(base_path)

    def galaxies(self, snapshot, fields=None):
        return self.catalog_reader.load_subhalos(snapshot, fields)

    def header_params(self, snapshot):
        return self.catalog_reader.load_header(snapshot)
