import astropy.units as u
import numpy as np
import h5py
import os
import eagleSqlTools
from constants import NFILES, DATA_PATH, SNAPSHOT_PATH, SNAPSHOT_PATH_SUFFIX


class EagleGalaxyDataReader(object):
    """
    Reads data for a specific galaxy from the EAGLE simulation.
    The galaxy is identified by the subhalo_id and the snapshot.
    """
    def __init__(self, gid, gn, sgn, centre, centre_vel, sim_reader):
        self.gid = gid
        self.gn = gn
        self.sgn = sgn
        self.centre = centre
        self.centre_vel = centre_vel

        self.simulation_reader = sim_reader
        self.simulation_reader.read_header()
        self.scaling_factor = self.simulation_reader.a
        self.hubble_param = self.simulation_reader.h
        self.omega0 = self.simulation_reader.omega0
        self.omegaLambda = self.simulation_reader.omegaLambda
        self.boxsize = self.simulation_reader.boxsize / self.simulation_reader.h * 1e3  # convert to kpc
        self.snapshot = self.simulation_reader.snapshot

        # Save the masks for this galaxy in a dictionary {itype: mask},
        # so that they are not re-calculated every time
        self.masks = {}

    def galaxy_mask(self, itype):
        #if itype in self.masks:
        #    return self.masks[itype]

        # Load data, then mask to selected GroupNumber and SubGroupNumber.
        gns = self.simulation_reader.read_dataset(itype, 'GroupNumber')
        sgns = self.simulation_reader.read_dataset(itype, 'SubGroupNumber')

        mask = np.logical_and(gns == self.gn, sgns == self.sgn)
        #self.masks[itype] = mask

        return mask

    def read_galaxy_coordinates(self, itype, wrap_around_centre=True):
        """ Resolve coordinates of itype particles for this galaxy.
        Coordinates are wrapped around the centre to account for periodicity, if required. """

        data = self.read_galaxy_attr(itype, 'Coordinates')
        data *= u.cm.to(u.kpc)

        # Periodic wrap coordinates around centre.
        if wrap_around_centre:
            boxsize = self.boxsize
            data = np.mod(data - self.centre + 0.5 * boxsize, boxsize) + self.centre - 0.5 * boxsize

        data = data - self.centre
        return data

    def read_galaxy_particle_ids(self, itype):
        """ Resolve particle ids of itype particles for this galaxy. """
        data = self.read_galaxy_attr(itype, 'ParticleIDs')
        return data

    def read_galaxy_velocities(self, itype, relative_to_centre=True):
        """ Resolve velocities of itype particles for this galaxy. """

        data = self.read_galaxy_attr(itype, 'Velocity')
        data *= u.cm.to(u.km)
        if relative_to_centre:
            data -= self.centre_vel
        return data

    def read_galaxy_mass(self, itype):
        """ Resolve mass of itype particles for this galaxy. """

        data = self.read_galaxy_attr(itype, 'Mass')
        data *= u.g.to(u.Msun)
        return data

    def read_galaxy_metallicity(self):
        """ Resolve metallicity of star particles for this galaxy. """

        data = self.read_galaxy_attr(4, 'Metallicity')
        return data

    def read_galaxy_stellar_formation_time(self):
        """ Resolve formation time of star particles for this galaxy
        (expansion factor when this star was born). """

        data = self.read_galaxy_attr(4, 'StellarFormationTime')
        return data

    def read_galaxy_gas_density(self):
        """ Resolve density of gas particles for this galaxy. """

        data = self.read_galaxy_attr(0, 'Density')
        return data

    def read_galaxy_gas_temperature(self):
        """ Resolve density of gas particles for this galaxy. """

        data = self.read_galaxy_attr(0, 'Temperature')
        return data

    def read_galaxy_sfr(self):
        """ Resolve star formation rate of each gas particle for this galaxy. """

        data = self.read_galaxy_attr(0, 'StarFormationRate')
        return data

    def read_galaxy_attr(self, itype, attr):
        """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        extract the requested attribute of all particles of a selected type. """

        # Load data, then mask to selected GroupNumber and SubGroupNumber.
        mask = self.galaxy_mask(itype)
        data = self.simulation_reader.read_dataset(itype, attr)[mask]

        return data


class EagleSnapshotDataReader(object):
    """
    Reads data for a specific snapshot for the EAGLE simulation.
    """

    def __init__(self, snapshot):
        path = SNAPSHOT_PATH % SNAPSHOT_PATH_SUFFIX[snapshot]
        self.path = path
        self.snapshot = snapshot

        # Initialize header variables
        self.a, self.h, self.boxsize = None, None, None
        self.omega0, self.omegaLambda, self.omegaBaryon = None, None, None

        self.dataset_dict = {}
        self.star_cutout_path = os.path.join(DATA_PATH, 'Maps_star_cutouts/'
                                                        'star_cutout_{}.h5py'.format(snapshot))

    def save_cutout(self):
        """ Save a cutout of star particles for speed-up """
        if os.path.exists(self.star_cutout_path):
            print('Star cutout %s already exists!' % self.star_cutout_path)
            return

        output_file = h5py.File(self.star_cutout_path, 'w')
        stars_group = output_file.create_group('STARS')

        fields = ('GroupNumber', 'SubGroupNumber', 'Coordinates', 'Mass',
                  'ParticleIDs', 'Velocity', 'Metallicity', 'StellarFormationTime')

        for field in fields:
            data = self.read_dataset(4, field, initial=True)
            stars_group.create_dataset(field, data=data)

        print('Finished writing star cutout ', self.star_cutout_path)
        output_file.close()

    def read_dataset(self, itype, att, nfiles=NFILES, keep_in_memory=False, initial=False):
        """ Read a selected dataset, itype is the PartType and att is the attribute name. """

        if itype == 4 and not initial:
            with h5py.File(self.star_cutout_path, 'r') as cutout_file:
                data = np.array(cutout_file['STARS/{}'.format(att)])
                return data

        if (itype, att) in self.dataset_dict:
            return self.dataset_dict[(itype, att)]

        # special case for DM mass
        if itype == 1 and att == 'Mass':
            data = self.read_dataset_dm_mass()
            if keep_in_memory:
                self.dataset_dict[(itype, att)] = data
            return data

        # Output array.
        data = []

        # Loop over each file and extract the data.
        for i in range(nfiles):
            f = h5py.File(self.path.format(i), 'r')

            tmp = f['PartType%i/%s' % (itype, att)][...]
            data.append(tmp)

            # Get conversion factors.
            cgs = f['PartType%i/%s' % (itype, att)].attrs.get('CGSConversionFactor')
            aexp = f['PartType%i/%s' % (itype, att)].attrs.get('aexp-scale-exponent')
            hexp = f['PartType%i/%s' % (itype, att)].attrs.get('h-scale-exponent')

            # Get expansion factor and Hubble parameter from the header.
            a = f['Header'].attrs.get('Time')
            h = f['Header'].attrs.get('HubbleParam')

            f.close()

        # Combine to a single array.
        if len(tmp.shape) > 1:
            data = np.vstack(data)
        else:
            data = np.concatenate(data)

        # Convert to physical.
        if data.dtype != np.int32 and data.dtype != np.int64:
            data = np.multiply(data, cgs * a ** aexp * h ** hexp, dtype='f8')

        if keep_in_memory:
            self.dataset_dict[(itype, att)] = data

        return data

    def read_header(self):
        """ Read various attributes from the header group. """
        f = h5py.File(self.path.format(0), 'r')
        self.a = f['Header'].attrs.get('Time')  # Scale factor.
        self.h = f['Header'].attrs.get('HubbleParam')  # h.
        self.boxsize = f['Header'].attrs.get('BoxSize')  # L [Mph/h].

        # total matter density in units of the critical density
        self.omega0 = f['Header'].attrs.get('Omega0')
        # density parameter corresponding to the cosmological constant
        self.omegaLambda = f['Header'].attrs.get('OmegaLambda')
        # mean baryon density in units of the critical density
        self.omegaBaryon = f['Header'].attrs.get('OmegaBaryon')
        f.close()

        return

    def read_dataset_dm_mass(self):
        """ Special case for the mass of dark matter particles. """
        f = h5py.File(self.path.format(0), 'r')
        h = f['Header'].attrs.get('HubbleParam')
        a = f['Header'].attrs.get('Time')
        dm_mass = f['Header'].attrs.get('MassTable')[1]
        n_particles = f['Header'].attrs.get('NumPart_Total')[1]

        # Create an array of length n_particles each set to dm_mass.
        m = np.ones(n_particles, dtype='f8') * dm_mass

        # Use the conversion factors from the mass entry in the gas particles.
        cgs = f['PartType0/Mass'].attrs.get('CGSConversionFactor')
        aexp = f['PartType0/Mass'].attrs.get('aexp-scale-exponent')
        hexp = f['PartType0/Mass'].attrs.get('h-scale-exponent')
        f.close()

        # Convert to physical.
        m = np.multiply(m, cgs * a ** aexp * h ** hexp, dtype='f8')

        return m


class EagleCatalogReader(object):
    """
    Reads catalog data for the EAGLE simulation using the eagleSqlTools module to
    connect to the database with your username and password.
    If the password is not given, the module will prompt for it.
    """
    def __init__(self, username, password):

        self.con = eagleSqlTools.connect(username, password=password)

    def query(self, sim_name, query, count=2, extra=()):
        # Execute query.
        try:
            args = count * (sim_name, ) + extra
            print(query,  args)
            query_data = eagleSqlTools.execute_query(self.con, query % args)
        except Exception as e:
            print(e)
            raise e

        return query_data


class EagleReader(object):
    """
    Reads data for EAGLE simulation.
    """
    def __init__(self, db_user, db_password, snapshot):
        self.catalog_reader = EagleCatalogReader(db_user, db_password)
        self.snapshot_reader = EagleSnapshotDataReader(snapshot)

    def query(self, sim_name, query, count=1, extra=()):
        return self.catalog_reader.query(sim_name, query, count=count, extra=extra)




