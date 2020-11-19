"""
 * Create RFT files for observation using input file and
 * complete RFT data using RMS data (well trajectories and grid)

Result:
  *.txt files written to ert/input/observations
  *.obs files for each well written to ert/input/observations
  well_date_rft.txt
"""
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline

CONFIG = {
        "rmsproject": project, # ROXAPI magic
        "gridname": "Simgrid",
        "zonename": "Zone",
        "interpolation": "cubic",
    "absolute_error": True,
    "obs_error": 3,
    "rftcsvfile": "../../ert/input/observations/rft/rft_input_table_4modelcase.csv",
    "exportdir": "../../ert/input/observations/rft",
    "aliasfile": "../input/well_modelling/well_info/rms_eclipse.csv",
    "aliascolumn_rms": "RMS_WELL_NAME",
    "aliascolumn_ecl": "ECLIPSE_WELL_ANME"

################################

# User settings
################################

grid_n = "Simgrid"
zone_n = "Zone"

# Interpolation method: linear if True, cubic spline if False
lin_interp = False

# Default observation error (if not specified). Given as an absolute value (3 bars by default).
# It could also be defined as a ratio of the observation value (2% or 5% for example).
# absolute_err: True for absolute error, False for relative error as percentage of the RFT value.
absolute_error = True
obs_error = 3.0

# RFT input file
# input_file = '../../ert/input/observations/rft/rft_input_table_4truth.csv'
input_file = "../../ert/input/observations/rft/rft_input_table_4modelcase.csv"
# CHANGE TO ISO8601 IN THIS FILE!!!!

# Location where to export the observation output files
path = "../../ert/input/observations/rft/"

# Well name conversion (RMS to Eclipse)
alias_file = "../input/well_modelling/well_info/rms_eclipse.csv"
rms_name = "RMS_WELL_NAME"  # Name of the column containing RMS well names
ecl_name = "ECLIPSE_WELL_NAME"  # Name of the column containing Eclipse well names

rft_prefix = ""  # Could be 'RFT_', 'R_', 'R', etc.
# To be used if the well_name in input_file correspond to the standard well name in RMS
# instead of the duplicated well associated with RFT data.
# If the name used in the input_file already includes it, then use: rft_prefix = ''
# The final RMS well_name will be considered to be "rft_prefix + <well_name from input_file>"


def make_alias_dict(alias_file, rms_name, ecl_name):
    """
    Create a correspondance dictionary so that well_dict[ <RMS wellname> ] = <Eclipse wellname>
    """
    df = pd.read_csv(alias_file, index_col=rms_name)
    well_dict = df.to_dict()
    well_dict = well_dict[ecl_name]

    return well_dict


def interp_from_md(wname, md_value, linear):
    """
    Function to interpolate East, North, TVD values of a well point
    defined by its name (wname) and its corresponding MD point (md_value).

    The interpolation of the well trajectory will be linear if linear = True,
    using cubic spline otherwise.
    """

    # Get trajectory survey points coordinates = [MD, EAST, NORTH, TVD]
    traj = project.wells[wname].wellbore.trajectories["Drilled trajectory"]
    coord = traj.survey_point_series.get_measured_depths_and_points()

    # Interpolate to get corresponding EAST, NORTH and TVD
    if linear:
        print("   Interpolating (linear) East, North, TVD from MD")
        # Values are interpolated using a linear interpolation
        itp_x = np.interp(md_value, coord[:, 0], coord[:, 1])  # East
        itp_y = np.interp(md_value, coord[:, 0], coord[:, 2])  # North
        itp_z = np.interp(md_value, coord[:, 0], coord[:, 3])  # TVD
    else:
        print("   Interpolating (cubic spline) East, North, TVD from MD")
        # Values are interpolated using CubicSpline
        cs_x = CubicSpline(coord[:, 0], coord[:, 1])  # East
        cs_y = CubicSpline(coord[:, 0], coord[:, 2])  # North
        cs_z = CubicSpline(coord[:, 0], coord[:, 3])  # TVD

        itp_x = float(cs_x(md_value))
        itp_y = float(cs_y(md_value))
        itp_z = float(cs_z(md_value))

    return itp_x, itp_y, itp_z


# def xyz_to_md(wellname, xyz, interpolation="cubic"):
def interp_from_xyz(wname, x0, y0, z0, linear):
    """Interpolate MD value of a well point defined by its name (wellname)
    and its corresponding East, North, TVD values (x0, y0, z0).

    The interpolation of the TVD well trajectory will be used
    (with linear algorithm if linear = True, using cubic spline otherwise)
    if the well is strictly going downward (TVD series stricly increasing).

    If the well is horizontal or in a hook shape, a projection of the point
    on the well trajectory will be used, using cubic spline interpolation
    of the well trajectory points.
    """

    # Get trajectory survey points coordinates = [MD, EAST, NORTH, TVD]
    traj = project.wells[wname].wellbore.trajectories["Drilled trajectory"]
    coord = traj.survey_point_series.get_measured_depths_and_points()

    # Interpolate to get corresponding MD

    print("   Interpolating MD from East, North, TVD")

    # Check if the well is strictly going down.
    strictly_downward = True
    for j in range(0, len(coord[:, 3]) - 1):
        if coord[j + 1, 3] < coord[j, 3]:
            strictly_downward = False

    if strictly_downward:

        # Interpolation from TVD survey can be used
        if linear:
            md_value = np.interp(z0, coord[:, 3], coord[:, 0])
        else:
            md_value = float(CubicSpline(coord[:, 3], coord[:, 0])(z0))
        print("   Methode 1: MD = {}".format(md_value))

    else:

        # In that case, the previous method is not valid,
        # MD must be used as a parametrization of
        # the trajectory coordinates (East, North, TVD)

        # RFT pressure point X0(x0,y0,z0) defined by:
        #     x0 = obs['EAST']
        #     y0 = obs['NORTH']
        #     z0 = obs['TVD']

        # distance of point X0 to any point X(x,y,z) of the well trajectory is:
        #     dist(X,X0)  = sqrt[ (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ]
        #     dist(MD,X0) = sqrt[ (f1(MD)-x0)**2 + (f2(MD)-y0)**2 + (f3(MD)-z0)**2 ]

        # We are looking for the point X, belonging to the well, which is the closest
        # from X0, so we want to find the value of MD for which minimizes dist(MD,X0).
        # This MD value is a root of the derivative of the distance. We drop the sqrt()
        # since this is equivalent to minimazing the square of the distance.
        #
        #     d dist(MD,X0) / d MD = 0
        #
        # <=> 2[ (f1(MD)-x0)*f1'(MD) + (f2(MD)-y0)*f2'(MD) + (f3(MD)-z0)*f3'(MD) ] = 0
        #
        # <=> (f1(MD)-x0)*f1'(MD) + (f2(MD)-y0)*f2'(MD) + (f3(MD)-z0)*f3'(MD) = 0
        #
        # using the notation f' for the derivative of function f
        #
        # This can be computed analytically when using polynomial interpolations for
        # f1, f2, f3 (having this function as a new polynom and finding its real roots),
        # but this method is quite unstable depending the degree of the polynom,
        # it may respect the survey points but being quite erratic in-between...
        #
        # Instead, we use cubic spline interpolations which provide a much smoother curve,
        # but we cannot define analytically the sum and multiplication of splines,
        # so we need to compute numerically the distance and find the minimum

        # x = f1(MD) for any point X(x,y,z) of the trajectory defined by MD (x = East)
        poly_x = CubicSpline(coord[:, 0], coord[:, 1])
        # y = f2(MD) for any point X(x,y,z) of the trajectory defined by MD (y = North)
        poly_y = CubicSpline(coord[:, 0], coord[:, 2])
        # z = f3(MD) for any point X(x,y,z) of the trajectory defined by MD (z = TVD)
        poly_z = CubicSpline(coord[:, 0], coord[:, 3])

        # compute the square of the distance every 5 cm MD (can take a couple of seconds)
        step = 0.05
        dist_min = np.inf
        closest_md = None
        for my_md in np.arange(0, coord[-1, 0], step):
            dist = (
                (poly_x(my_md) - x0) ** 2
                + (poly_y(my_md) - y0) ** 2
                + (poly_z(my_md) - z0) ** 2
            )
            if dist < dist_min:
                closest_md = my_md
                dist_min = dist
        dist_min = dist_min ** (0.5)
        md_value = round(closest_md, 2)
        print("   Methode 2: MD = {} (mismatch = {:.4f} m)".format(md_value, dist_min))

        # Alternative using polynoms
        # degree = 7
        # poly_x = np.polyfit(coord[:, 0], coord[:, 1], degree) # x = f1(MD)
        # poly_y = np.polyfit(coord[:, 0], coord[:, 2], degree) # y = f2(MD)
        # poly_z = np.polyfit(coord[:, 0], coord[:, 3], degree) # z = f3(MD)

        # pder_x = np.polyder(poly_x)             # derivative, f1'(MD)
        # pder_y = np.polyder(poly_y)             # derivative, f2'(MD)
        # pder_z = np.polyder(poly_z)             # derivative, f3'(MD)

        # poly_x[-1] = poly_x[-1] - x0            # f1(MD) := f1(MD)-x0
        # poly_y[-1] = poly_y[-1] - y0            # f2(MD) := f2(MD)-y0
        # poly_z[-1] = poly_z[-1] - z0            # f3(MD) := f3(MD)-z0

        # Derivative of the distance function
        # pder_dist  = np.polyadd(np.polyadd(np.polymul(poly_x, pder_x),
        #                                    np.polymul(poly_y, pder_y)),
        #                         np.polymul(poly_z, pder_z))

        # Minimum of the (squared) distance function is the real root of this polynom
        # roots = np.roots(pder_dist)
        # interp_md = roots[np.isreal(roots)].real

        # Only one root should be real
        # assert len(interp_md) == 1, ('0 or more than 1 point have been found ({}) to fit '
        #                              'the well trajectory for this RFT data point!'
        #                              .format(len(interp_md)))

        # md_value = interp_md[0]
        # print('   Methode 3: MD = {}'.format(md_value))

    return md_value


###############################################################################
#                                                                             #
#                                 MAIN SCRIPT                                 #
#                                                                             #
###############################################################################

data = pd.read_csv(input_file)

well_names = data["WELL_NAME"].unique().tolist()
grid_model = project.grid_models[grid_n]
grid = grid_model.get_grid()

alias_list = make_alias_dict(alias_file, rms_name, ecl_name)

datfilename = path + "well_date_rft.txt"
datfile = open(datfilename, "w")

for well in well_names:

    well_input = well
    well = rft_prefix + well

    assert well in project.wells.keys(), 'No well "{}"!'.format(well)

    welldata = data[(data.WELL_NAME == well_input)]
    welldata.reset_index(drop=False)
    # number of observations defined for the given well
    n_obs = welldata.shape[0]

    print("Well {}: {} observation(s)".format(well, n_obs))
    obs_list = []
    for i in range(0, n_obs):
        obs = welldata.iloc[i].to_dict()

        is_md = not np.isnan(obs["MD"])
        tvd = not np.isnan(obs["TVD"])
        if np.isnan(obs["EAST"]) or np.isnan(obs["NORTH"]) or np.isnan(obs["TVD"]):
            xyz = False
        else:
            xyz = True

        if not xyz:
            assert is_md, (
                "MD must be specified as input " "if (East, North, TVD) is not given!"
            )
            # Compute X,Y,Z
            itpx, itpy, itpz = interp_from_md(well, obs["MD"], lin_interp)
            obs["EAST"] = itpx
            obs["NORTH"] = itpy
            obs["TVD"] = itpz

        elif not is_md:
            # Compute MD
            obs["MD"] = interp_from_xyz(
                well, obs["EAST"], obs["NORTH"], obs["TVD"], lin_interp
            )

        cell_index = grid.get_cells_at_points([obs["EAST"], obs["NORTH"], obs["TVD"]])
        assert not (cell_index is None), (
            "The RFT point is outside of the grid! "
            "(MD = {:.2f}, East = {:.1f}, North = {:.1f}, "
            "TVD = {:.1f})".format(obs["MD"], obs["EAST"], obs["NORTH"], obs["TVD"])
        )

        # If a zone is provided:
        # check if the cell is in the right zone, if not: warning
        if not pd.isnull(obs["ZONE"]):
            # Find zone corresponding to the EAST, NORTH, TVD

            zone_val = grid_model.properties[zone_n].get_values()[cell_index]
            zone_str = grid_model.properties[zone_n].code_names[zone_val]

            if obs["ZONE"] == zone_str:
                print("   Zone check: OK")
            else:
                print(
                    "   WARNING: input zone ({}) does not match the "
                    "grid zone ({}) at given coordinates!".format(obs["ZONE"], zone_str)
                )

        assert not np.isnan(obs["VALUE"]), (
            "The observation value must " "be specified as input!"
        )

        # If no error measurement is provided, provide default error based on inputs
        if np.isnan(obs["ERROR"]):
            if absolute_error:
                # if absolute value
                obs["ERROR"] = obs_error
            else:
                # if relative ratio
                obs["ERROR"] = obs["VALUE"] * obs_error

        obs_list.append(obs)

    # Generate output text files from the different obs dict() for each well
    if well in alias_list:
        print("rms well name:", well, "|  eclipse well name:", alias_list[well])
        well = alias_list[well]
    else:
        print("------eclipse alias name was not found for", well)
        print("------please update the alias file", alias_file)

    txtfilename = path + well + ".txt"
    obsfilename = path + well + ".obs"

    obsfile = open(obsfilename, "w")
    txtfile = open(txtfilename, "w")

    for obs in obs_list:
        obsfile.write("{:.3f}   {:.3f}\n".format(obs["VALUE"], obs["ERROR"]))
        txtfile.write(
            "{:.3f}   {:.3f}   {:.3f}   {:.3f}   {}\n".format(
                obs["EAST"], obs["NORTH"], obs["MD"], obs["TVD"], obs["ZONE"]
            )
        )

    print("   Data written to " + obsfilename)
    print("   Data written to " + txtfilename)

    obsfile.close()
    txtfile.close()

    # TODO: DATE should be ISO-8601 in csv file, convert to ERT's DD.MM.YYYY here.
    unique_dates = welldata["DATE"].unique().tolist()
    for count, date in enumerate(unique_dates, 1):
        date_with_spaces = date.replace(".", " ")
        datfile.write("{}   {}   {}\n".format(well, date_with_spaces, count))

datfile.close()
print("Data written to " + datfilename)

print("Done.")
