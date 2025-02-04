import hyram.phys.api as phys_api
import numpy as np
from hyram.phys import Orifice, Flame, Jet

release_fluid = phys_api.create_fluid('H2',
                                        temp=288,  # K
                                        pres=35e6,  # Pa
                                        phase='none')

ambient_fluid = phys_api.create_fluid('AIR',temp=288,  # K
                                        pres=101325)  # Pa
leak_diam = 0.003  # m
orifice = Orifice(leak_diam)
flame = Flame(release_fluid,
                orifice,
                ambient_fluid,
                verbose=True)

result = phys_api.analyze_jet_plume(ambient_fluid, release_fluid, leak_diam, mass_flow=None,
                      rel_angle=(0*np.pi/2), dis_coeff=1, nozzle_model='yuce',
                      create_plot=True, contours=None,
                      xmin=None, xmax=None, ymin=None, ymax=None, vmin=0, vmax=0.1,
                      plot_title="Mole Fraction of Leak",
                      filename=None, output_dir=None, verbose=False)

# Select initial conditions of relase 
# Select composition of release
# Select ambient conditions
# Select leak size and angle
# Select nozzle model
# For jet plume plot the concentration contours
# Calculate unconfined over pressure 

