import math

import numpy as np

import pytest

from fmu.tools.rms import create_rft_ertobs


class SurveyPointSeries:
    def __init__(self):
        pass

    def get_measured_depths_and_points(self):
        return np.array([[0, 0, 0, 0], [1, 0, 0, 1]])


class Trajectory:
    def __init__(self):
        self.survey_point_series = SurveyPointSeries()


class Wellbore:
    def __init__(self):
        self.trajectories = {"Drilled trajectory": Trajectory()}


class Well:
    def __init__(self):
        self.wellbore = Wellbore()


class Grid:
    def __init__(self):
        pass

    def get_cells_at_points(xyz):
        return 13345


class Property:
    def __init__(self):
        self.code_names = {"Valyzar": 1}

    def get_values(self):
        pass


class GridModel:
    def __init__(self):
        self.properties = {"Valyzar": Property()}

    def get_grid(self):
        return Grid()


class RMSMockedProject:
    def __init__(self):
        self.wells = {
            "R_A2": Well(),
            "R_A3": Well(),
            "R_A4": Well(),
            "R_A5": Well(),
            "R_A6": Well(),
        }
        self.grid_models = {"Simgrid": GridModel()}


def test_get_well_coords():
    """Reliance on this test depends on whether the mocked
    RMS project resembles the real ROXAPI.
    """
    rms_project_mock = RMSMockedProject()

    assert (
        create_rft_ertobs.get_well_coords(rms_project_mock, "R_A4")
        == np.array([[0, 0, 0, 0], [1, 0, 0, 1]])
    ).all()


@pytest.mark.parametrize(
    "coords, expected",
    [
        (np.array([[0, 0, 0, 0], [1, 0, 0, 1]]), True),
        (np.array([[0, 0, 0, 0], [1, 0, 0, 1], [2, 0, 0, 0]]), False),
        (np.array([[0, 0, 0, 0], [1, 0, 0, 1], [1.001, 0, 0, 0.999]]), False),
    ],
)
def test_strictly_downward(coords, expected):
    assert create_rft_ertobs.strictly_downward(coords) == expected


def test_interp_from_md():
    coords = np.array([[0, 0, 0, 0], [1, 0, 0, 1]])
    assert create_rft_ertobs.interp_from_md(0.5, coords, interpolation="linear") == (
        0,
        0,
        0.5,
    )
    assert create_rft_ertobs.interp_from_md(0.5, coords, interpolation="cubic") == (
        0,
        0,
        0.5,
    )


def test_interp_from_xyz():
    coords = np.array([[0, 0, 0, 0], [1, 0, 0, 1]])
    assert create_rft_ertobs.interp_from_xyz((0, 0, 0.5), coords) == 0.5

    # A wellpath going straight down, and then up again at 45 degrees:
    coords = np.array([[0, 0, 0, 0], [1, 0, 0, 1], [1 + math.sqrt(2), 1, 0, 0]])
    assert (
        create_rft_ertobs.interp_from_xyz((0, 0, 0.5), coords, interpolation="linear")
        == 0.35  # verify this number
    )
    assert (
        create_rft_ertobs.interp_from_xyz(
            (0.5, 0.0, 0.5), coords, interpolation="linear"
        )
        == 2.0  # verify this number
    )
