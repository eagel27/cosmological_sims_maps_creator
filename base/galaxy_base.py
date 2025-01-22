import os.path
import sys
import h5py
import pandas as pd
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib
import matplotlib.ticker as mticker
from collections import Counter
from scipy.ndimage import gaussian_filter
from multiprocessing import Pool
from .constants import PLANE_LOS_PROJECTION, LOS_PROJECTION, SOLAR_METALLICITY
from .alignments import ALIGNMENTS
from .utils import calc_percentiles


class GalaxyBase(object):

    def __init__(self, should_load_datasets=True):
        self.star_coordinates = None
        self.star_particle_ids = None
        self.star_mass = None
        self.star_velocities = None
        self.star_metallicity = None
        self.star_formation_time = None
        self.vel_system = None
        self.spin = None
        self.cosm_params = {}
        self.galaxy_reader = None

        self.view_id = 0
        self.alignments = []

        if should_load_datasets:
            self.load_datasets()

    def get_id(self):
        pass

    def load_datasets(self):
        self.star_coordinates = self.galaxy_coordinates(4)
        self.star_particle_ids = self.galaxy_particle_ids(4)
        self.star_mass = self.galaxy_mass(4)
        self.star_velocities = self.galaxy_velocities(4)
        self.star_metallicity = self.galaxy_metallicity()
        self.star_formation_time = self.galaxy_stellar_formation_time()
        self.vel_system = self.calc_systemic_velocity()
        self.spin = self.calculate_spin()
        self.cosm_params = self.get_cosm_params()

    def get_cosm_params(self):
        pass

    def get_exsitu_df(self, snapshot):
        pass

    def galaxy_mass(self, itype):
        return self.galaxy_reader.read_galaxy_mass(itype)

    def galaxy_particle_ids(self, itype):
        return self.galaxy_reader.read_galaxy_particle_ids(itype)

    def galaxy_velocities(self, itype):
        return self.galaxy_reader.read_galaxy_velocities(itype)

    def galaxy_coordinates(self, itype):
        return self.galaxy_reader.read_galaxy_coordinates(itype)

    def galaxy_metallicity(self):
        return self.galaxy_reader.read_galaxy_metallicity()

    def galaxy_stellar_formation_time(self):
        return self.galaxy_reader.read_galaxy_stellar_formation_time()

    def galaxy_sfr(self):
        return self.galaxy_reader.read_galaxy_sfr()

    def galaxy_gas_density(self):
        return self.galaxy_reader.read_galaxy_gas_density()

    def galaxy_gas_temperature(self):
        return self.galaxy_reader.read_galaxy_gas_temperature()

    def calculate_age(self, expansionFactor):
        hubble_param = self.cosm_params.get('hubble_param')
        omegaLambda = self.cosm_params.get('omegaLambda')
        omega0 = self.cosm_params.get('omega0')

        hubble_time = 9.778 / hubble_param  # Gyr

        factor1 = 2.0 / 3.0 / np.sqrt(omegaLambda)
        term1 = np.sqrt(omegaLambda / omega0) * expansionFactor ** 1.5
        term2 = np.sqrt(1 + omegaLambda / omega0 * expansionFactor ** 3)
        factor2 = np.log(term1 + term2)

        time = factor1 * factor2 * hubble_time

        return time

    def calc_systemic_velocity(self):
        pass

    def calculate_spin(self):
        coords = self.star_coordinates
        masses = self.star_mass
        velocities = self.star_velocities

        velocities = velocities - self.vel_system

        central_mask = np.linalg.norm(coords, axis=1) <= 30
        spin_numeric = np.average(np.cross(coords[central_mask],
                                           velocities[central_mask]),
                                  weights=masses[central_mask], axis=0)
        return spin_numeric

    def edge_on_rotator(self, star_parameter):
        """
        Rotates the coordinates or velocities of the stars according to the galactic spin provided.
        The spin of the galaxy is aligned with the Z axis of the simulation.

        Parameters:
            -spin_gal: 1 x 3 array
                Spin of the galaxy in the simulation rest frame

            -star_parameter:     N x 3 array
                Coordinates or velocities of the stars in the simulation rest frame

        Returns:
            -Star_parameter_rotated:  N x 3 array
                Star parameter rotated according to the spin of the galaxy
        """
        spin_gal = self.spin
        return self.rotate_galaxy(spin_gal, star_parameter)

    def rotate_x_axis(self, theta):
        f1 = [1, 0, 0]
        f2 = [0, np.cos(theta), -np.sin(theta)]
        f3 = [0, np.sin(theta), np.cos(theta)]
        rotation_matrix = np.array([f1, f2, f3])
        return rotation_matrix

    def rotate_y_axis(self, theta):
        f1 = [np.cos(theta), 0, np.sin(theta)]
        f2 = [0, 1, 0]
        f3 = [-np.sin(theta), 0, np.cos(theta)]
        rotation_matrix = np.array([f1, f2, f3])
        return rotation_matrix

    def rotate_z_axis(self, theta):
        f1 = [np.cos(theta), -np.sin(theta), 0]
        f2 = [np.sin(theta), np.cos(theta), 0]
        f3 = [0, 0, 1]
        rotation_matrix = np.array([f1, f2, f3])
        return rotation_matrix

    def rotate_galaxy(self, rot_vector, star_parameter):
        """
        Rotates the coordinates or velocities of the stars according to rotation vector provided.
        The rotation vector is aligned with the Z axis of the simulation.
        """
        spin_x_angle = np.arctan2(rot_vector[1], rot_vector[2])
        x_rotation_matrix = self.rotate_x_axis(spin_x_angle)
        spin_x_rotated = np.dot(x_rotation_matrix, rot_vector)

        spin_y_angle = np.arctan2(spin_x_rotated[0], spin_x_rotated[2])

        y_rotation_matrix = self.rotate_y_axis(-spin_y_angle)
        total_rotation_matrix = np.dot(y_rotation_matrix, x_rotation_matrix)
        star_parameter_rotated = np.matmul(total_rotation_matrix, star_parameter.transpose()).transpose()
        return star_parameter_rotated

    def change_galaxy_view(self, star_parameter):
        if not self.alignments:
            return star_parameter
        view = np.array(self.alignments[self.view_id])
        dist = 2 * self.cosm_params['boxsize']
        rot_vector = view * dist + self.centre
        return self.rotate_galaxy(rot_vector, star_parameter)

    def binned_mapping(self, plane_xy, values, statistic='sum'):
        stamp_side_pix = 129  # pixels of map per side
        stamp_rel_side = 4  # n_times the half stellar mass radius (hmsr) of the galaxy
        hmsr_kpc = self.hmrg

        bins = np.linspace(- hmsr_kpc * stamp_rel_side / 2,
                           hmsr_kpc * stamp_rel_side / 2, stamp_side_pix)

        hist = np.histogram2d(plane_xy[:, 1], plane_xy[:, 0], bins=[bins, bins])

        mmap, _, _, _ = stats.binned_statistic_2d(plane_xy[:, 1], plane_xy[:, 0], values,
                                                  statistic=statistic,
                                                  bins=[bins, bins])

        a = Counter(hist[0].flat)
        percentiles = calc_percentiles(a)

        lim = max(percentiles[75][1], 2)
        mmap[np.where(hist[0] <= lim)] = np.nan

        return mmap

    def project_coordinates(self):
        coords = self.star_coordinates

        if self.align:
            coords = self.edge_on_rotator(coords)
        else:
            coords = self.change_galaxy_view(coords)

        coords_projected = coords[:, PLANE_LOS_PROJECTION[self.los]]
        return coords_projected

    def project_velocities(self, vel_central_calc=True):
        velocities = self.star_velocities

        # Use centre of galaxy as the rest frame
        velocities = velocities - self.vel_system

        if self.align:
            velocities = self.edge_on_rotator(velocities)
        else:
            velocities = self.change_galaxy_view(velocities)

        velocities_projected = velocities[:, LOS_PROJECTION[self.los]]

        return velocities_projected

    def stellar_mass_map(self):
        coords_projected = self.project_coordinates()
        masses = self.star_mass
        mmap = self.binned_mapping(coords_projected, masses)
        return mmap

    def stellar_velocity_map(self):
        coords_projected = self.project_coordinates()
        velocities_projected = self.project_velocities()
        mmap = self.binned_mapping(coords_projected, velocities_projected, statistic='mean')
        return mmap

    def stellar_velocity_dispersion_map(self):
        coords_projected = self.project_coordinates()
        velocities_projected = self.project_velocities()
        mmap = self.binned_mapping(coords_projected, velocities_projected, statistic='std')
        return mmap

    def stellar_metallicity_map(self, mass_map):
        coords_projected = self.project_coordinates()
        metallicity = self.star_metallicity / SOLAR_METALLICITY
        masses = self.star_mass
        weighted_metallicity = np.multiply(masses, metallicity)
        mmap = self.binned_mapping(coords_projected, weighted_metallicity)
        mmap = np.divide(mmap, mass_map)
        return mmap

    def stellar_formation_time_map(self, mass_map):
        coords_projected = self.project_coordinates()
        formation_time = self.star_formation_time
        masses = self.star_mass
        weighted_formation_time = np.multiply(masses, formation_time)
        mmap = self.binned_mapping(coords_projected, weighted_formation_time)
        mmap = np.divide(mmap, mass_map)
        return mmap

    def stellar_age_map(self, mass_map):
        coords_projected = self.project_coordinates()
        formation_time = self.star_formation_time
        stellar_ages = self.calculate_age(self.cosm_params['scaling_factor']) - self.calculate_age(formation_time)
        masses = self.star_mass
        weighted_stellar_ages = np.multiply(masses, stellar_ages)
        mmap = self.binned_mapping(coords_projected, weighted_stellar_ages)
        mmap = np.divide(mmap, mass_map)
        return mmap

    def create_2d_maps(self, alignments):
        views_maps = []
        for view_id in range(len(alignments) or 1):
            self.view_id = view_id
            self.alignments = alignments
            mass_map = self.stellar_mass_map()
            vel_map = self.stellar_velocity_map()
            vel_disp_map = self.stellar_velocity_dispersion_map()
            metal_map = self.stellar_metallicity_map(mass_map)
            sft_map = self.stellar_formation_time_map(mass_map)
            age_map = self.stellar_age_map(mass_map)
            views_maps.append((mass_map, vel_map, vel_disp_map, metal_map, sft_map, age_map))
        return self.get_id(), views_maps

    def save_2d_maps(self, alignments):
        galaxy_id, (mass_map, vel_map, vel_disp_map,
                    metal_map, sft_map, age_map) = self.create_2d_maps(alignments)

        galaxy_group = self.output_file.create_group('{id}_data'.format(id=int(galaxy_id)))

        galaxy_maps_group = galaxy_group.create_group('MAPS')

        galaxy_maps_group.create_dataset('mass_map', data=mass_map)
        galaxy_maps_group.create_dataset('vel_map', data=vel_map)
        galaxy_maps_group.create_dataset('vel_disp_map', data=vel_disp_map)
        galaxy_maps_group.create_dataset('metal_map', data=metal_map)
        galaxy_maps_group.create_dataset('sft_map', data=sft_map)
        galaxy_maps_group.create_dataset('age_map', data=age_map)

    def plot_2d_maps(self, maps_file, subhalo_mass, subhalo_exsitu, snapshot, path, view,
                     fig=None, subfigs=None):
        plt.style.use('science')
        # plt.rcParams.update({'font.size': 20})
        # plt.rc('axes', titlesize=20, labelsize=20)

        path = os.path.join(path, 'EAGLE_TNGlike')
        if not os.path.exists(path):
            os.makedirs(path)

        galaxy_id = '{}_{}'.format(self.get_id(), snapshot)
        mass_map = maps_file['{}_{}_data/MAPS/mass_map'.format(galaxy_id, view)][...]
        vel_map = maps_file['{}_{}_data/MAPS/vel_map'.format(galaxy_id, view)][...]
        vel_disp_map = maps_file['{}_{}_data/MAPS/vel_disp_map'.format(galaxy_id, view)][...]
        metal_map = maps_file['{}_{}_data/MAPS/metal_map'.format(galaxy_id, view)][...]
        sft_map = maps_file['{}_{}_data/MAPS/sft_map'.format(galaxy_id, view)][...]
        age_map = maps_file['{}_{}_data/MAPS/age_map'.format(galaxy_id, view)][...]

        sigma = 0.5
        i = 1
        if fig is None:
            i = 0
            fig = plt.figure(figsize=(50, 20), dpi=100)
            subfigs = fig.subfigures(nrows=2, ncols=1)

        subfig = subfigs[i]
        ax = subfig.add_subplot(2, 5, 1 + 5 * i)
        mass_map = np.log10(mass_map)
        mass_map[~np.isfinite(mass_map)] = 0
        mass_map = gaussian_filter(mass_map, sigma=sigma)

        c = ax.imshow(mass_map, 'jet', rasterized=False, origin='lower', vmin=0, vmax=12)
        plt.colorbar(c)
        ax.set_title('Log(Mass)')

        ax = subfig.add_subplot(2, 5, 2 + 5 * i)
        vel_map[~np.isfinite(vel_map)] = 0
        vel_map = gaussian_filter(vel_map, sigma=sigma)
        lim = min(min(abs(np.nanmin(vel_map)), abs(np.nanmax(vel_map))), 200)
        vel_map = vel_map.clip(-lim, lim)
        c = ax.imshow(vel_map, 'jet', rasterized=False, origin='lower', vmin=-200, vmax=200)
        plt.colorbar(c)
        ax.set_title('Velocity (km/s)')

        ax = subfig.add_subplot(2, 5, 3 + 5 * i)
        vel_disp_map[~np.isfinite(vel_disp_map)] = 0
        vel_disp_map = gaussian_filter(vel_disp_map, sigma=sigma)
        c = ax.imshow(vel_disp_map, 'jet', rasterized=True, origin='lower', vmin=0, vmax=350)
        plt.colorbar(c)
        ax.set_title('Velocity dispersion (km/s)')

        ax = subfig.add_subplot(2, 5, 4 + 5 * i)
        metal_map = np.log10(metal_map)
        metal_map[~np.isfinite(metal_map)] = 0
        metal_map = gaussian_filter(metal_map, sigma=sigma)
        metal_map = metal_map.clip(-1, 2)
        c = ax.imshow(metal_map, 'jet', rasterized=True, origin='lower', vmin=-1, vmax=2)
        plt.colorbar(c)
        ax.set_title('Log(Metallicity)')

        # ax = fig.add_subplot(1, 5, 5)
        # sft_map[~np.isfinite(sft_map)] = 0
        # sft_map = gaussian_filter(sft_map, sigma=sigma)
        # c = ax.imshow(sft_map, 'jet', rasterized=True, origin='lower', interpolation='bilinear')
        # plt.colorbar(c)
        # ax.set_title('Stellar Formation Time')

        ax = subfig.add_subplot(2, 5, 5 + 5 * i)
        age_map[~np.isfinite(age_map)] = 0
        age_map = gaussian_filter(age_map, sigma=sigma)
        c = ax.imshow(age_map, 'jet', rasterized=True, origin='lower', interpolation='bilinear', vmin=0, vmax=11)
        plt.colorbar(c)
        ax.set_title('Stellar Age (Gyr)')

        # Create plot title
        f = mticker.ScalarFormatter(useOffset=False, useMathText=True)

        # g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.2e' % x))
        # fmt = mticker.FuncFormatter(g)
        def as_si(x, ndp):
            s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
            m, e = s.split('e')
            return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))
        plt.style.use('science')
        plt.tight_layout()
        plt.show()

    def get_central_sigma(self, maps_file, subhalo_exsitu, snapshot):
        galaxy_id = '{}_{}'.format(self.get_id(), snapshot)
        vel_disp_map = maps_file[galaxy_id + '_0_data/MAPS/vel_disp_map'][...]
        central_sigma = np.mean(vel_disp_map[59:69, 59:69])
        mass = subhalo_exsitu.TotalMass.iloc[0]
        return mass, central_sigma


class GalaxyMapsBase(object):
    MAPS_RESULTS_PATH = None  # overwrite
    SNAPSHOT = None  # overwrite
    SIM_NAME = None  # overwrite
    SIM = None  # overwrite
    CATALOG_DS_PATH = None  # overwrite

    def __init__(self, align=False, los='Z'):
        self.suffix = None

    def get_output_path(self):
        output_file_name = 'maps_{}_balanced{}.h5py'.format(
            self.SIM_NAME, '_{}'.format(self.suffix) if self.suffix else '')
        output_path = os.path.join(self.MAPS_RESULTS_PATH, output_file_name)
        if not os.path.exists(self.MAPS_RESULTS_PATH):
            os.makedirs(self.MAPS_RESULTS_PATH)
        return output_path

    def load_balanced_data(self):
        balanced_dataset = pd.read_hdf(self.CATALOG_DS_PATH, '/data')
        snapshots = balanced_dataset.Snapshot.unique()
        galaxies_data = []
        for snapshot in snapshots:
            snapshot_df = balanced_dataset[balanced_dataset.Snapshot == snapshot]
            snapshot_galaxy_ids = set(snapshot_df.GalaxyID)
            snapshot_catalog_data = self.load_catalog_data(snapshot)
            # Keep only data in balanced dataset
            for d in snapshot_catalog_data:
                if d[0] not in snapshot_galaxy_ids:
                    continue
                g_data = snapshot_df[snapshot_df.GalaxyID == d[0]].to_dict(orient='records')[0]
                g_data.update({'GalaxyData': list(d)})
                galaxies_data.append(g_data)

        return galaxies_data[self.id_start:self.id_stop]

    def load_catalog_data(self, snapshot):
        """
        :rtype: list containing the galaxies catalog data
        """
        return

    def create_galaxy_object(self, subhalo_data, snapshot, output_file=None, should_load_datasets=True):
        """
        :rtype: Galaxy object
        """
        return

    def create_galaxy_maps(self, subhalo_data_dict):
        subhalo_data = subhalo_data_dict['GalaxyData']
        snapshot = subhalo_data_dict['Snapshot']
        num_alignments = subhalo_data_dict['Num_of_Aligns']
        alignments = ALIGNMENTS[num_alignments]
        galaxy_object = self.create_galaxy_object(subhalo_data, snapshot)
        g_id, maps = galaxy_object.create_2d_maps(alignments)
        del galaxy_object
        return '{}_{}'.format(g_id, snapshot), maps

    def create_maps(self):
        galaxy_catalog_data = self.load_balanced_data()
        galaxies_count = len(galaxy_catalog_data)

        pool = Pool(processes=8)
        results = list(tqdm.tqdm(pool.imap(self.create_galaxy_maps, galaxy_catalog_data),
                                 total=galaxies_count))
        pool.close()
        pool.join()

        output_file = h5py.File(self.get_output_path(), 'w')
        for galaxy_id, view_maps in results:
            for view in range(len(view_maps)):
                maps = view_maps[view]
                galaxy_group = output_file.create_group('{gid}_{view}_data'.format(gid=galaxy_id,
                                                                                   view=view))
                (mass_map, vel_map, vel_disp_map,
                 metal_map, sft_map, age_map) = maps
                galaxy_maps_group = galaxy_group.create_group('MAPS')

                galaxy_maps_group.create_dataset('mass_map', data=mass_map)
                galaxy_maps_group.create_dataset('vel_map', data=vel_map)
                galaxy_maps_group.create_dataset('vel_disp_map', data=vel_disp_map)
                galaxy_maps_group.create_dataset('metal_map', data=metal_map)
                galaxy_maps_group.create_dataset('sft_map', data=sft_map)
                galaxy_maps_group.create_dataset('age_map', data=age_map)

        print('Finished writing to ', self.get_output_path())
        output_file.close()

    def save_2d_maps(self):
        galaxies_catalog_data = self.load_balanced_data()
        galaxies_count = len(galaxies_catalog_data)
        output_file = h5py.File(self.get_output_path(), 'r')

        path = os.path.join(self.MAPS_RESULTS_PATH, 'Plots')
        for i, subhalo_data_dict in enumerate(galaxies_catalog_data):
            print('Plotting for Galaxy {}/{}'.format(i + 1, galaxies_count))
            subhalo_data = subhalo_data_dict['GalaxyData']
            snapshot = subhalo_data_dict['Snapshot']
            num_alignments = subhalo_data_dict['Num_of_Aligns']
            subhalo_exsitu = subhalo_data_dict['ExSitu_Fraction']
            subhalo_mass = subhalo_data_dict['Mass']
            g_object = self.create_galaxy_object(subhalo_data,
                                                 snapshot,
                                                 output_file=output_file,
                                                 should_load_datasets=False)
            for view in range(num_alignments):
                g_object.plot_2d_maps(output_file, subhalo_mass, subhalo_exsitu, snapshot, path, view)

            del g_object

        output_file.close()
