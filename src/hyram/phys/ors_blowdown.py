import numpy as np
from scipy.constants import liter, bar, atm, psi, milli
from scipy import interpolate
from matplotlib import pyplot as plt
from hyram.phys import Source, Fluid, Orifice
from hyram.utilities import misc_utils


# Set up parameters
temp = 273  # K...assumed initial tank temp is ambient air temp. Table 1 of 'Exp. investigation of hydrogen release' by Ekoto
pressure = 232e5 + atm  # Pa...Table 1 of above source
orifice_diameter = 0.00356  # m...Table 1 of above source
discharge_coeff = 0.75  # ...Table 1 of above source
tank_volume = 0.363  # m^3...Table 1 of above source
rel_species = "CO2"

fluid = Fluid(species=rel_species, P=pressure, T=temp)

# plt.plot(calc_time, calc_flowrate)
# plt.plot(calc_time, source.pres)
# plt.xlabel("Time (s)")
# plt.ylabel("Flow Rate (kg/s)")
# plt.show()

hole_id = [2, 7, 36, 112]
vol = [0.5, 1, 10, 100]
results = []
figures = []
for i, v in enumerate(vol):
    plt.figure(i + 1)
    for j, hole in enumerate(hole_id):
        print(f"Calculating for {v} m3 and {hole} mm hole")
        orifice_diameter = hole * milli
        tank_volume = v
        orifice = Orifice(orifice_diameter, discharge_coeff)
        source = Source(tank_volume, fluid)
        calc_flowrate, _, calc_time, _ = source.empty(orifice, t_empty=900)
        if calc_time[-1] < 899:
            calc_time += [899, 900]
            calc_flowrate += [0, 0]
        # fig = source.plot_time_to_empty()
        results.append(source)
        plt.semilogy(calc_time, calc_flowrate, label=f"Hole size {hole} mm")
    plt.title(f"Volume {v} m3")
    plt.legend(loc="best")
    plt.xlabel("Time (s)")
    plt.ylabel("Release rate (kg/s)")
    plt.tight_layout()
    plt.savefig(f"Volume_{v}_m3.png", dpi=300)
    plt.show()
    # plt.plot(calc_time, calc_flowrate)
    # plt.plot(source.ts, source.pres)
    # plt.xlabel("Time (s)")
    # plt.ylabel("Flow Rate (kg/s)")
    # plt.show()
