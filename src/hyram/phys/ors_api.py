import win32api
import win32con
import xlwings as xw
import hyram.phys.api as phys_api
import numpy as np
from hyram.phys import Orifice, Flame, Jet, Source


def run_confined_overpressure(sheet_name, filename=None):
    if filename == None:
        wb = xw.Book.caller()
    else:
        wb = xw.Book(filename)

    wb.app.calculation = "manual"
    sheet = wb.sheets[sheet_name]

    species = sheet.range("REL_COMP").value

    if species == "mixture":
        species = {}
        for component, molefrac in zip(
            sheet.range("HYRAM_COMPS").value, sheet.range("HYRAM_MOLEFRACS").value
        ):
            if molefrac and molefrac > 0:
                species[component] = molefrac

        # molefracs =

    nozzle_model = sheet.range("NOZZLE_MODEL").value

    if nozzle_model == "YuceilOtugen":
        nozzle_model = "yuce"
    if nozzle_model == "EwanMoodie":
        nozzle_model = "ewan"
    if nozzle_model == "Birch":
        nozzle_model = "birc"
    if nozzle_model == "Birch2":
        nozzle_model = "bir2"
    if nozzle_model == "Molkov":
        nozzle_model = "molk"

    tank_volume = sheet.range("TANK_VOLUME").value
    rel_height = sheet.range("REL_HEIGHT").value
    enclos_height = sheet.range("ENCLOSURE_HEIGHT").value
    floor_ceil_area = sheet.range("FLOOR_CEIL_AREA").value
    floor_vent_xarea = sheet.range("FLOOR_VENT_AREA").value
    floor_vent_height = sheet.range("FLOOR_VENT_HEIGHT").value
    ceil_vent_xarea = sheet.range("CEIL_VENT_AREA").value
    ceil_vent_height = sheet.range("CEIL_VENT_HEIGHT").value
    max_time = sheet.range("MAX_TIME").value
    dist_rel_to_wall = sheet.range("DIST_REL_TO_WALL").value

    if sheet.range("BLOWDOWN").value == "YES" and tank_volume > 0:
        is_steady = False
    else:
        is_steady = True
    # win32api.MessageBox(wb.app.hwnd, str(sheet.range("BLOWDOWN").value), 'Warning', win32con.MB_ICONINFORMATION)

    times = np.linspace(0, max_time, int(max_time))

    release_fluid = phys_api.create_fluid(
        species,
        temp=sheet.range("STAG_TEMP").value + 273.15,  # K
        pres=sheet.range("STAG_PRES").value * 1e5,  # Pa
        phase="none",
    )

    ambient_fluid = phys_api.create_fluid(
        "AIR",
        temp=sheet.range("AMB_TEMP").value + 273.15,  # K
        pres=sheet.range("AMB_PRES").value * 1e5,
    )  # Pa
    leak_diam = sheet.range("LEAK_DIAMETER").value / 1000  # m
    rel_angle = sheet.range("RELEASE_ANGLE").value / 180 * np.pi
    Cd = sheet.range("LEAK_CD").value

    result = phys_api.analyze_accumulation(
        ambient_fluid,
        release_fluid,
        tank_volume,
        leak_diam,
        rel_height,
        enclos_height,
        floor_ceil_area,
        ceil_vent_xarea,
        ceil_vent_height,
        floor_vent_xarea,
        floor_vent_height,
        times,
        orif_dis_coeff=Cd,
        ceil_vent_coeff=1,
        floor_vent_coeff=1,
        vol_flow_rate=0,
        dist_rel_to_wall=dist_rel_to_wall,
        tmax=max_time,
        rel_area=None,
        rel_angle=rel_angle,
        nozzle_key=nozzle_model,
        x0=0,
        y0=0,
        nmax=1000,
        temp_pres_points=None,
        pres_ticks=None,
        is_steady=is_steady,
        create_plots=True,
        output_dir="out",
        verbose=False,
    )

    sheet.pictures.add(
        result["pres_plot_filepath"],
        name="Overpressure",
        update=True,
        left=sheet.range("C66").left,
        top=sheet.range("C66").top,
        scale=0.9,
    )
    sheet.pictures.add(
        result["layer_plot_filepath"],
        name="Layer",
        update=True,
        left=sheet.range("C93").left,
        top=sheet.range("C93").top,
        scale=0.9,
    )

    sheet.pictures.add(
        result["trajectory_plot_filepath"],
        name="Trajectory",
        update=True,
        left=sheet.range("C121").left,
        top=sheet.range("C121").top,
        scale=0.9,
    )
    sheet.pictures.add(
        result["mass_plot_filepath"],
        name="Mass",
        update=True,
        left=sheet.range("C145").left,
        top=sheet.range("C145").top,
        scale=0.9,
    )
    # sheet.pictures.add(result['mass_flow_plot_filepath'], name='Blowdown', update=True,
    #                        left=sheet.range('B122').left, top=sheet.range('B122').top, scale=1.3)
    if not is_steady:
        sheet.pictures.add(
            result["mass_flow_plot_filepath"],
            name="Blowdown",
            update=True,
            left=sheet.range("B172").left,
            top=sheet.range("B172").top,
            scale=1.3,
        )
    else:
        try:
            sheet.pictures("Blowdown").delete()
        except:
            pass

    sheet.range("MAX_OVERPRESSURE").value = result["overpressure"] / 1e5


def run_blowdown(sheet_name, filename=None):
    if filename == None:
        wb = xw.Book.caller()
    else:
        wb = xw.Book(filename)

    wb.app.calculation = "manual"
    sheet = wb.sheets[sheet_name]

    species = sheet.range("REL_COMP").value

    if species == "mixture":
        species = {}
        for component, molefrac in zip(
            sheet.range("HYRAM_COMPS").value, sheet.range("HYRAM_MOLEFRACS").value
        ):
            if molefrac and molefrac > 0:
                species[component] = molefrac

        # molefracs =

    tank_volume = sheet.range("TANK_VOLUME").value
    max_time = sheet.range("MAX_TIME").value

    win32api.MessageBox(
        wb.app.hwnd, str(max_time), "Warning", win32con.MB_ICONINFORMATION
    )

    release_fluid = phys_api.create_fluid(
        species,
        temp=sheet.range("STAG_TEMP").value + 273.15,  # K
        pres=sheet.range("STAG_PRES").value * 1e5,  # Pa
        phase="none",
    )

    ambient_fluid = phys_api.create_fluid(
        "AIR",
        temp=sheet.range("AMB_TEMP").value + 273.15,  # K
        pres=sheet.range("AMB_PRES").value * 1e5,
    )  # Pa
    leak_diam = sheet.range("LEAK_DIAMETER").value / 1000  # m
    Cd = sheet.range("LEAK_CD").value

    orifice = Orifice(leak_diam, Cd)
    source = Source(tank_volume, release_fluid)
    calc_flowrate, _, calc_time, _ = source.empty(orifice, t_empty=max_time)

    fig = source.plot_time_to_empty()

    sheet.pictures.add(
        fig,
        name="Blowdown",
        update=True,
        left=sheet.range("B68").left,
        top=sheet.range("B68").top,
        scale=1.3,
    )

    sheet.range("MAX_RELEASE").value = max(calc_flowrate)
    sheet.range("MIN_RELEASE").value = min(calc_flowrate)
    sheet.range("DURATION").value = max(calc_time)

    times = np.linspace(0, int(max(calc_time)), 20)
    pres = np.interp(times, source.ts, source.pres)
    mdot = np.interp(times, source.ts, source.mdot)
    mass = np.interp(times, source.ts, source.sol[0])

    sheet.range("C33").options(transpose=True).value = times
    sheet.range("D33").options(transpose=True).value = pres / 1e5
    sheet.range("E33").options(transpose=True).value = mdot
    sheet.range("F33").options(transpose=True).value = mass


def run_unconfined_overpressure(sheet_name, filename=None):
    if filename == None:
        wb = xw.Book.caller()
    else:
        wb = xw.Book(filename)

    wb.app.calculation = "manual"
    sheet = wb.sheets[sheet_name]

    species = sheet.range("REL_COMP").value

    if species == "mixture":
        species = {}
        for component, molefrac in zip(
            sheet.range("HYRAM_COMPS").value, sheet.range("HYRAM_MOLEFRACS").value
        ):
            if molefrac and molefrac > 0:
                species[component] = molefrac

        # molefracs =

    nozzle_model = sheet.range("NOZZLE_MODEL").value

    calculation_method = sheet.range("CALCULATION_METHOD").value
    # locations = [(float(sheet.range("OPERP_X").value), float(sheet.range("OPERP_Y").value), float(sheet.range("OPERP_Z").value))]
    coord_x = sheet.range("OVERP_X").value
    coord_y = sheet.range("OVERP_Y").value
    coord_z = sheet.range("OVERP_Z").value
    locations = [(coord_x, coord_y, coord_z)]

    bst_flame_speed = 0.35
    if calculation_method == "BST":
        bst_flame_speed = sheet.range("BST_FLAME_SPEED").value

    if nozzle_model == "YuceilOtugen":
        nozzle_model = "yuce"
    if nozzle_model == "EwanMoodie":
        nozzle_model = "ewan"
    if nozzle_model == "Birch":
        nozzle_model = "birc"
    if nozzle_model == "Birch2":
        nozzle_model = "bir2"
    if nozzle_model == "Molkov":
        nozzle_model = "molk"

    release_fluid = phys_api.create_fluid(
        species,
        temp=sheet.range("STAG_TEMP").value + 273.15,  # K
        pres=sheet.range("STAG_PRES").value * 1e5,  # Pa
        phase="none",
    )

    ambient_fluid = phys_api.create_fluid(
        "AIR",
        temp=sheet.range("AMB_TEMP").value + 273.15,  # K
        pres=sheet.range("AMB_PRES").value * 1e5,
    )  # Pa
    leak_diam = sheet.range("LEAK_DIAMETER").value / 1000  # m
    rel_angle = sheet.range("RELEASE_ANGLE").value / 180 * np.pi
    Cd = sheet.range("LEAK_CD").value

    result = phys_api.compute_overpressure(
        calculation_method,
        locations,
        ambient_fluid,
        release_fluid,
        leak_diam,
        mass_flow=None,
        release_angle=rel_angle,
        discharge_coefficient=Cd,
        nozzle_model=nozzle_model,
        bst_flame_speed=bst_flame_speed,
        tnt_factor=0.03,
        flammability_limits=None,
        origin_at_orifice=False,
        create_overpressure_plot=True,
        create_impulse_plot=True,
        overpressure_plot_filename=None,
        impulse_plot_filename=None,
        output_dir="out",
        verbose=False,
        overp_contours=None,
        overp_xlims=None,
        overp_ylims=None,
        overp_zlims=None,
        impulse_contours=None,
        impulse_xlims=None,
        impulse_ylims=None,
        impulse_zlims=None,
    )

    sheet.pictures.add(
        result["overp_plot_filepath"],
        name="Overpressure",
        update=True,
        left=sheet.range("B32").left,
        top=sheet.range("B32").top,
    )

    if calculation_method != "Bauwens":
        sheet.pictures.add(
            result["impulse_plot_filepath"],
            name="Impulse",
            update=True,
            left=sheet.range("B69").left,
            top=sheet.range("B69").top,
        )
    else:
        sheet.pictures("Impulse").delete()

    sheet.range("OVERP_MASS_FLOW").value = result["mass_flow_rate"]
    sheet.range("OVERP_PRESSURE").value = result["overpressures"] / 1000
    sheet.range("OVERP_IMPULSE").value = result["impulses"] / 1000
    sheet.range("OVERP_MASS").value = result["flam_or_det_mass"]


def run_jet_flame(sheet_name, filename=None):
    if filename == None:
        wb = xw.Book.caller()
    else:
        wb = xw.Book(filename)

    wb.app.calculation = "manual"
    sheet = wb.sheets[sheet_name]

    species = sheet.range("REL_COMP").value

    if species == "mixture":
        species = {}
        for component, molefrac in zip(
            sheet.range("HYRAM_COMPS").value, sheet.range("HYRAM_MOLEFRACS").value
        ):
            if molefrac and molefrac > 0:
                species[component] = molefrac

        # molefracs =

    nozzle_model = sheet.range("NOZZLE_MODEL").value
    rel_hum = sheet.range("AMB_HUM").value
    wind_speed = sheet.range("WIND_SPEED").value

    if nozzle_model == "YuceilOtugen":
        nozzle_model = "yuce"
    if nozzle_model == "EwanMoodie":
        nozzle_model = "ewan"
    if nozzle_model == "Birch":
        nozzle_model = "birc"
    if nozzle_model == "Birch2":
        nozzle_model = "bir2"
    if nozzle_model == "Molkov":
        nozzle_model = "molk"

    release_fluid = phys_api.create_fluid(
        species,
        temp=sheet.range("STAG_TEMP").value + 273.15,  # K
        pres=sheet.range("STAG_PRES").value * 1e5,  # Pa
        phase="none",
    )

    ambient_fluid = phys_api.create_fluid(
        "AIR",
        temp=sheet.range("AMB_TEMP").value + 273.15,  # K
        pres=sheet.range("AMB_PRES").value * 1e5,
    )  # Pa
    leak_diam = sheet.range("LEAK_DIAMETER").value / 1000  # m
    rel_angle = sheet.range("RELEASE_ANGLE").value / 180 * np.pi
    Cd = sheet.range("LEAK_CD").value

    coord_x = sheet.range("RADIATION_X").value
    coord_y = sheet.range("RADIATION_Y").value
    coord_z = sheet.range("RADIATION_Z").value

    result = phys_api.jet_flame_analysis(
        ambient_fluid,
        release_fluid,
        leak_diam,
        mass_flow=None,
        dis_coeff=Cd,
        rel_angle=rel_angle,
        nozzle_key=nozzle_model,
        rel_humid=rel_hum,
        create_temp_plot=True,
        temp_plot_filename=None,
        temp_plot_title="",
        temp_contours=None,
        temp_xlims=None,
        temp_ylims=None,
        analyze_flux=True,
        flux_plot_filename=None,
        flux_coordinates=[(coord_x, coord_y, coord_z)],
        flux_contours=None,
        flux_xlims=None,
        flux_ylims=None,
        flux_zlims=None,
        output_dir="out",
        verbose=False,
        wind_speed=wind_speed,
    )

    sheet.pictures.add(
        result[1],
        name="Flame heat flux",
        update=True,
        left=sheet.range("B32").left,
        top=sheet.range("B32").top,
    )

    sheet.pictures.add(
        result[0],
        name="Flame temperature",
        update=True,
        left=sheet.range("B68").left,
        top=sheet.range("B68").top,
    )

    sheet.range("FLAME_MASS_FLOW").value = result[3]
    sheet.range("POSITIONAL_FLUX").value = result[2]
    sheet.range("RADIATIVE_POWER").value = result[4] / 1000
    sheet.range("VISIBLE_LENGTH").value = result[5]
    sheet.range("RADIANT_FRAC").value = result[6]


def run_jet_plume(sheet_name, filename=None):
    if filename == None:
        wb = xw.Book.caller()
    else:
        wb = xw.Book(filename)

    wb.app.calculation = "manual"
    sheet = wb.sheets[sheet_name]

    species = sheet.range("REL_COMP").value

    if species == "mixture":
        species = {}
        for component, molefrac in zip(
            sheet.range("HYRAM_COMPS").value, sheet.range("HYRAM_MOLEFRACS").value
        ):
            if molefrac and molefrac > 0:
                species[component] = molefrac

        # molefracs =

    nozzle_model = sheet.range("NOZZLE_MODEL").value

    if sheet.range("CONTOUR").value:
        contour = sheet.range("CONTOUR").value
    else:
        contour = None
    if sheet.range("XMIN").value is not None:
        xmin = sheet.range("XMIN").value
    else:
        xmin = None
    if sheet.range("XMAX").value is not None:
        xmax = sheet.range("XMAX").value
    else:
        xmax = None
    if sheet.range("YMIN").value is not None:
        ymin = sheet.range("YMIN").value
    else:
        ymin = None
    if sheet.range("YMAX").value is not None:
        ymax = sheet.range("YMAX").value
    else:
        ymax = None

    if sheet.range("VMAX").value:
        vmax = sheet.range("VMAX").value
        vmin = 0
    else:
        vmax = 0.1
        vmin = 0

    if nozzle_model == "YuceilOtugen":
        nozzle_model = "yuce"
    if nozzle_model == "EwanMoodie":
        nozzle_model = "ewan"
    if nozzle_model == "Birch":
        nozzle_model = "birc"
    if nozzle_model == "Birch2":
        nozzle_model = "bir2"
    if nozzle_model == "Molkov":
        nozzle_model = "molk"

    release_fluid = phys_api.create_fluid(
        species,
        temp=sheet.range("STAG_TEMP").value + 273.15,  # K
        pres=sheet.range("STAG_PRES").value * 1e5,  # Pa
        phase="none",
    )

    ambient_fluid = phys_api.create_fluid(
        "AIR",
        temp=sheet.range("AMB_TEMP").value + 273.15,  # K
        pres=sheet.range("AMB_PRES").value * 1e5,
    )  # Pa
    leak_diam = sheet.range("LEAK_DIAMETER").value / 1000  # m
    rel_angle = sheet.range("RELEASE_ANGLE").value / 180 * np.pi
    Cd = sheet.range("LEAK_CD").value

    orifice = Orifice(leak_diam)
    # flame = Flame(release_fluid, orifice, ambient_fluid, verbose=True)

    result = phys_api.analyze_jet_plume(
        ambient_fluid,
        release_fluid,
        leak_diam,
        mass_flow=None,
        rel_angle=rel_angle,
        dis_coeff=Cd,
        nozzle_model=nozzle_model,
        create_plot=True,
        contours=contour,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
        # ymin=-10,
        # ymax=0,
        vmin=vmin,
        vmax=vmax,
        plot_title="Mole Fraction of Leak",
        filename=None,
        output_dir=None,
        verbose=False,
    )
    # sheet.pictures['Jet Plume'].delete
    sheet.pictures.add(
        result["plot"],
        name="Jet Plume",
        update=True,
        left=sheet.range("B32").left,
        top=sheet.range("B32").top,
    )

    sheet.range("RES_MASS_FLOW").value = result["mass_flow_rate"]

    LFL = [key for key in result["mole_frac_dists"].keys()][0]
    sheet.range("RES_MIN_X").value = result["mole_frac_dists"][LFL][0][0]
    sheet.range("RES_MAX_X").value = result["mole_frac_dists"][LFL][0][1]
    sheet.range("RES_MIN_Y").value = result["mole_frac_dists"][LFL][1][0]
    sheet.range("RES_MAX_Y").value = result["mole_frac_dists"][LFL][1][1]
    sheet.range("LFL").value = LFL


# Select initial conditions of relase
# Select composition of release
# Select ambient conditions
# Select leak size and angle
# Select nozzle model
# For jet plume plot the concentration contours
# Calculate unconfined over pressure

if __name__ == "__main__":
    release_fluid = phys_api.create_fluid(
        "CO2", temp=288, pres=10e6, phase="none"  # K  # Pa
    )

    ambient_fluid = phys_api.create_fluid("AIR", temp=288, pres=101325)  # K  # Pa
    leak_diam = 0.01  # m
    orifice = Orifice(leak_diam)
    # flame = Flame(release_fluid, orifice, ambient_fluid, verbose=True)

    result = phys_api.analyze_jet_plume(
        ambient_fluid,
        release_fluid,
        leak_diam,
        mass_flow=None,
        rel_angle=(1 * np.pi / 2),
        dis_coeff=1,
        nozzle_model="yuce",
        create_plot=True,
        contours=[0.05],
        xmin=-2,
        xmax=2,
        ymin=0,
        ymax=15,
        vmin=0,
        vmax=0.1,
        plot_title="Mole Fraction of Leak",
        filename=None,
        output_dir=None,
        verbose=False,
    )
    # result = phys_api.jet_flame_analysis(ambient_fluid, release_fluid, leak_diam, mass_flow=None, dis_coeff=1,
    #    rel_angle=0,
    #    nozzle_key='yuce', rel_humid=0.89,
    #    create_temp_plot=True, temp_plot_filename=None, temp_plot_title="", temp_contours=None,
    #    temp_xlims=None, temp_ylims=None,
    #    analyze_flux=True, flux_plot_filename=None, flux_coordinates= [(0, 0, 0), (1, 1, 0)], flux_contours=None,
    #    flux_xlims=None, flux_ylims=None, flux_zlims=None,
    #    output_dir=None, verbose=False)
