import IDtrack

import Person

class Location:


    def __init__(self,x,y):
        self.location = (x,y)
        self.occupants = []


    def simulateSecond(self):

        for p in self.occupants:
            pass