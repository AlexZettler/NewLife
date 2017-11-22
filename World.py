import IDtrack
import Location
import Person
import sys
import random
import SomeStatStuff



class World:


    def __init__(self, savePath, generate=False):

        self.locations = IDtrack.IDtracker()
        self.sentientBeings = IDtrack.IDtracker()

        self.lookup = LookupTable(lookup)



        if generate:
            self.GenerateWorld(savePath,3)

        else:
            self.LoadWorld(savePath)

    def simulateSecond(self):

        for l in self.locations:

            l.simulateSecond()

    def WorldExists(self, savePath):
        pass

    def GenerateWorld(self,savePath,nodes):

        for n in range(nodes):

            npos = SomeStatStuff.genRandInsideUnitCircle(10)

            self.locations.store_val(npos)



    def LoadWorld(self,savePath):
        pass


class LookupTable:

    def __init__(self):

        self.RaceLookup = Person.RaceTree()
        self.ActionLookup = Person.ActionTree()


