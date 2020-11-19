class Wellbore:
    def __init__(self):
        self.trajectories = {"Drilled trajectory": "foo"}


class Well:
    def __init__(self):
        self.wellbore = Wellbore()



class RMSMockedProject:
    def __init__(self):
        self.wells = {"A-4": Well()}

project = RMSMockedProject()

project.wells["A-4"].wellbore.trajectories["Drilled trajectory"]
