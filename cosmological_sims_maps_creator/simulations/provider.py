from cosmological_sims_maps_creator.simulations.EAGLE import EagleGalaxyMaps
from cosmological_sims_maps_creator.simulations.EAGLE_IllustrisTNG import EagleIllustrisTNGGalaxyMaps
from cosmological_sims_maps_creator.simulations.IllustrisTNG import IllustrisTNGGalaxyMaps


class GalaxySimulationProvider:
    def __init__(
        self, galaxy_type, mode, align, id_start=None, id_stop=None, suffix=""
    ):
        self.galaxy_type = galaxy_type
        assert self.galaxy_type in ["EAGLE", "EAGLE_IllustrisTNG", "IllustrisTNG"]

        self.mode = mode
        self.align = align
        self.id_start = id_start
        self.id_stop = id_stop
        self.suffix = suffix

    def _get_galaxy_maps(self):
        if self.galaxy_type == "EAGLE":
            return EagleGalaxyMaps
        elif self.galaxy_type == "EAGLE_IllustrisTNG":
            return EagleIllustrisTNGGalaxyMaps
        elif self.galaxy_type == "IllustrisTNG":
            return IllustrisTNGGalaxyMaps

    def run(self):
        if self.align == 1:
            galaxy_maps = self._get_galaxy_maps()(
                align=True,
                los="X",
                suffix=self.suffix,
                id_start=self.id_start,
                id_stop=self.id_stop,
            )
        else:
            galaxy_maps = self._get_galaxy_maps()(
                suffix=self.suffix, id_start=self.id_start, id_stop=self.id_stop
            )

        if self.mode == 1:
            galaxy_maps.create_maps()
        elif self.mode == 2:
            galaxy_maps.save_2d_maps()
