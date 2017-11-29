import random
import xml.etree.ElementTree

import XML_Parse_Helper as ph


#


class Person:
    def __init__(self, name, race, interests, relations):

        self.name = name
        self.race = race

        self.keyAbilities = None

        self.interests = interests
        self.relations = relations

        self.happiness = 100
        self.hunger = 0

        self.timeSinceRest = 0
        self.rested = 100

        self.inventory = []



    def __str__(self):
        return " ".join(self.name)

    @classmethod
    def NewPerson(self):

        race = Person.RandomRace({"Human": 1})
        name = Person.GenName(race)
        interests = None
        relations = None

        return Person(name,race,interests,relations)


    @classmethod
    def RandomRace(cls, raceWeights):

        cumulativeWeight = 0
        tempRaceWeights = raceWeights

        for w in raceWeights.keys():

            #pulls all text attributes aka race names from the race lookup
            if w in [r.get("name") for r in RaceLookup.races.findall("Race")]:

                cumulativeWeight += raceWeights[w]

                #This could be a problem if the race inputs are not valid <POTENTIAL ERROR>
                tempRaceWeights[w] = cumulativeWeight

            else:
                raise KeyError

        #print(cumulativeWeight)
        randint = random.randint(1, cumulativeWeight)
        #print(randint)

        for r in tempRaceWeights.keys():
            if randint <= tempRaceWeights[r]:
                return r


        #Brute force search therefore error if it reaches here
        raise LookupError



    def GenerateInterests(self,race):
        raise NotImplementedError

    @classmethod
    def GenName(self, Race):

        #Foreach Race in xml file
        allRaces = RaceLookup.races.findall("Race")
        for r in allRaces:

            #Get element's text content
            if Race == r.get("name"):
                try:
                    nameList = r.find("Names")

                    firstNames = nameList.find("First").findall("li")
                    #print("There are "+str(len(firstNames))+" first names found for race "+Race)
                    firstName = firstNames[random.randint(0, len(firstNames))-1]
                    #print(firstName.text)

                    lastNames = nameList.find("Last").findall("li")
                    #print("There are "+str(len(lastNames))+" last names found for race "+Race)
                    lastName = lastNames[random.randint(0, len(lastNames))-1]
                    #print(lastName.text)

                    return (firstName.text, lastName.text)

                except:
                    raise LookupError

class ActionController:

    class Action:

        def __init__(self,name,target):
            try:
                pass
                #ActionTree.actions.findall("Action")

    def __init__(self):

        #action priorities as self.ActionPriorities[actionname] = weight
        self.ActionPriorities = {}

        self.currentAction = None
        self.currentActionTime = None

    def SimulateSecond(self):
        if self.currentActionTime >= ActionTree.actions.findall(self.currentAction)


class RaceTree:
    def __init__(self):
        self.races = xml.etree.ElementTree.parse('Races.xml').getroot()

        r = self.races.findall("Race")
        for sr in r:
            #Get element's text content
            print(sr.get("name"))

    @classmethod
    def raceValid(self,race):

        for r in RaceLookup.races.findall("Race"):

            #Get element's text content
            if race == r.get("name"):
                return True
        return False


class Interest:
    def __init__(self):

        self.interests = xml.etree.ElementTree.parse('Interests.xml').getroot()
        r = self.interests.findall("Interest")

class Desire:
    def __init__(self):
        self.desires = xml.etree.ElementTree.parse('Desire.xml').getroot()
        r = self.desires.findall("Desire")

class ActionTree:
    def __init__(self):
        self.actions = xml.etree.ElementTree.parse('Actions.xml').getroot()
        r = self.actions.findall("Action")


def test():
    print()
    print("*"*20)
    print("Test Starting")
    #print("*"*20)


    #print(Person.NewRandom({"Human":5,"Dwarf":1,"Elf":1,"Ork":1,"Minotaur":1,"Goblin":1}))

    #print("*"*20)
    print("Test Complete")
    print("*"*20)



#test()





#Look this over in the morning you dipshit... and go mow the grass when you fuck shit up some more



p = []
for i in range(10):
    p.append(Person.NewPerson())
    print(p[i])

#


