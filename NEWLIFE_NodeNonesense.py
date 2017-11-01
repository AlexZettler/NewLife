import math
import random
import copy

import IDtrack

import matplotlib.pyplot as pplot
import numpy

class NodeWrangler:
    def __init__(self):

        self.nodes = IDtrack.IDtracker()
        self.links = IDtrack.IDtracker()

        self._originalNodes = None

    def __str__(self):
        return "\n".join(str(self.nodes[n]) for n in self.nodes)


    def genNodes(self,n,radius=10):
        for i in range(n):
            rn = genRandInsideUnitCircle(radius)
            self.nodes.StoreVal(PhysicalNode(rn, 1.0, None))
        self._originalNodes = copy.deepcopy(self.nodes)
    def genAllPhysicalLinks(self,strength):
        #nkt = 0
        #print("Node key length: {}".format(len(self.nodeKeys)))
        #print(self.nodes.keys())

        completedkeys = []

        for key in self.nodes:

            completedkeys.append(key)
            connectToList = {k: self.nodes[k] for k in set(self.nodes.keys()) - set(completedkeys)}

            for connection in connectToList:

                pl = physicalNodeLink(self.nodes[key],self.nodes[connection])

                self.links.StoreVal(pl)

                self.nodes[key].links.append(pl)
                self.nodes[connection].links.append(pl)


    @classmethod
    def angFromSinCos(self,sinang,cosang):

        # Quadrant breakup
        #Q1
        if sinang > 0 and cosang >= 0:
            return cosang
        # Q2
        elif sinang < 0 and cosang >= 0:
            return math.pi - cosang
        # Q3
        elif sinang < 0 and cosang <= 0:
            return math.pi + cosang
        # Q4
        elif sinang < 0 and cosang >= 0:
            return 2 * math.pi + cosang

        else:
            raise ValueError

    def resetMoveVectors(self):
        for n in self.nodes:
            #print("Node: {}".format(self.nodes[n].id))
            #print("movVec: {}".format(type(self.nodes[n].moveVec)))
            self.nodes[n].moveVec = Pos2D((0, 0))

    def iterativeCalculateForceVectors(self,iterations,startStr,endStr):

        for i in range(iterations):
            self.CalculateForceVectors((i+1)/iterations * startStr/endStr)
            print("i/iterations: {}, startStr/endStr: {}, mul={}".format(i/iterations, startStr/endStr, i/iterations* startStr/endStr))
            self.ApplyVectors()


    def CalculateForceVectors(self, strength):
        print("Iterating through vectors for force calculation")

        for l in self.links:

            #(l.n1 * l.n2)/l.lengthSquared()

            hypSq = self.links[l].lengthSquared()
            hyp = math.sqrt(hypSq)

            sinang = math.sin(self.links[l].deltaY / hyp)
            cosang = math.cos(self.links[l].deltaX / hyp)
            ang = NodeWrangler.angFromSinCos(sinang,cosang)

            #print("cosang: {} sinang: {}".format(cosang*180/math.pi,sinang*180/math.pi))
            #n.moveVec.x += (l.n1 * l.n2)


            #need to recalculate hypsquared here so optional parameter to ease calculation
            chargemul = self.links[l].chargedParticleForce(hyp) * strength

            #print("N1:{},N2:{} ".format(self.links[l].n1.moveVec,self.links[l].n2.moveVec))

            self.links[l].n1.moveVec.x += chargemul / self.links[l].n1.weight * math.acos(cosang)
            self.links[l].n1.moveVec.y += chargemul / self.links[l].n1.weight * math.asin(sinang)

            self.links[l].n2.moveVec.x += -chargemul / self.links[l].n2.weight * math.acos(cosang)
            self.links[l].n2.moveVec.y += -chargemul / self.links[l].n2.weight * math.asin(sinang)



    def ApplyVectors(self):

        #Apply the moveVectors
        for n in self.nodes:

            if self.nodes[n].rooted == False:

                #print("type: {}, node x: {}, y: {}".format(type(self.nodes[n]),self.nodes[n].x, self.nodes[n].y)) #+= self.nodes[n].moveVec.x
                #print("type: {}, movevec x: {}, y: {}".format(type(self.nodes[n].moveVec),self.nodes[n].moveVec.x, self.nodes[n].moveVec.y))#+= self.nodes[n].moveVec.y
                self.nodes[n].x += self.nodes[n].moveVec.x
                self.nodes[n].y += self.nodes[n].moveVec.y

        self.resetMoveVectors()



    def dispPlot(self):

        #n = numpy.
        #frame = pplot.scatter(200,200)


        #numpy.random.seed(19680801)
        #data = numpy.random.randn(2, 100)

        xdata = [self.nodes[n].x for n in self.nodes]
        ydata = [self.nodes[n].y for n in self.nodes]

        xdataOrig = [self._originalNodes[n].x for n in self._originalNodes]
        ydataOrig = [self._originalNodes[n].y for n in self._originalNodes]


        data = numpy.array([xdata,ydata])


        #configures a grid of plots
        #fig, axs = pplot.subplots(2, 2, figsize=(5, 5)

        #

        #axs[0, 0].hist(data[0])
        pplot.scatter(data[0], data[1])
        #axs[0, 1].plot(data[0], data[1])
        #axs[1, 1].hist2d(data[0], data[1])

        pplot.show()


    def dispStreamLineGraph(self):


        xdata = [self.nodes[n].x for n in self.nodes]
        ydata = [self.nodes[n].y for n in self.nodes]

        xmoveData = [self.nodes[n].moveVec.x for n in self.nodes]
        ymoveData = [self.nodes[n].moveVec.y for n in self.nodes]

        data = numpy.array([xdata, ydata,xmoveData,ymoveData])


        #X, Y, U, V = zip(data)


        print("some graph data:\n\t{},{} pieces of data in xdata,ydata\n\t{},{} pieces of data in xmovedata,ymovedata".format(len(xdata),len(ydata),len(xmoveData),len(ymoveData)))

        #posData = numpy.array([xdata,ydata])
        #velData = numpy.array([xmoveData, ymoveData])

        pplot.quiver(data[0], data[1], data[2], data[3])
        pplot.show()

class Pos2D:

    def __init__(self,pos):
        #self.pos = None
        self.__pos = None
        self.pos = pos

            #print("Instantiating Pos2D - self.pos type:{}".format(str(type(self.pos))))

    def __str__(self):
        return "Pos2D-x:{}, y:{}".format(self.x, self.y)

    def __getitem__(self, item):
        if isinstance(item,int):
            try:
                return self.__pos[item]
            except:
                raise StopIteration

    def __setitem__(self, key, value):
        if isinstance(key,int) and key in (0, 1):
            self.__pos[key] = value

    @property
    def pos(self):
        return self.__pos

    @pos.setter
    def pos(self,val):
       #print("Pos set as:{}".format(val))

        if len(val) == 2 and type(val) in (list, set, tuple, Pos2D):
            self.__pos = [val[0], val[1]]
        else:
            raise ValueError


    @property
    def x(self):
        return self.__pos[0]

    @x.setter
    def x(self, val):

        try:
            #print("setting self.pos from {} to {}".format(self.pos[0],val))
            self.__pos[0] = val
        except:
            #print("Error setting self.pos of type:value from {}:{} to {}:{} - self.pos of type {}".format(type(self.pos[0]),self.pos[0], type(val),val,type(self.pos)))
            raise TypeError

    @property
    def y(self):
        return self.__pos[1]

    @y.setter
    def y(self, val):
        try:
            #print("setting self.pos from {} to {}".format(self.pos[0],val))
            self.__pos[1] = val
        except:
            #print("Error setting self.pos of type:value from {}:{} to {}:{} - self.pos of type {}".format(type(self.pos[0]),self.pos[0], type(val),val,type(self.pos)))
            raise TypeError

class Node(Pos2D):

    def __init__(self,pos,id):

        super().__init__(pos)

        #ID gets set from IDtracker
        self.id = id
        self.links = []

    def __str__(self):
        if type(id) != None:
            return "id: {}, x: {}, y: {}".format(self.id,self.x,self.y)
        else:
            return str(super())

class PhysicalNode(Node):

    def __init__(self, pos, weight, id):
        super().__init__(pos,id)

        self.weight = weight
        self.moveVec = Pos2D((0, 0))
        self.rooted = False

    def __str(self):
        return "{}, moveVec: {}".format(str(super), self.moveVec)


class nodeLink():

    def __init__(self,n1,n2):
        self.n1, self.n2 = n1, n2

    def lengthSquared(self):
        return self.deltaX**2 + self.deltaY**2

    def length(self):
        return math.sqrt(self.lengthSquared())

    def __len__(self):
        return self.length()

    @property
    def deltaX(self):
        return self.n2.x - self.n1.x

    @property
    def deltaY(self):
        return self.n2.y - self.n1.y




class physicalNodeLink(nodeLink):
    def __init__(self,n1, n2, force=1.0):
        super().__init__(n1, n2)

    def chargedParticleForce(self, lensquared):
        if type(lensquared)==float:
            return (self.n1.weight * self.n2.weight) / lensquared
        else:
            return ValueError

#@classmethod
def genRandInsideUnitCircle(radius=1.0):
    x = (random.random() - 0.5) * radius
    y = (random.random() - 0.5) * radius

    while x ** 2 + y ** 2 > radius ** 2:
        x = (random.random()-0.5) * radius
        y = (random.random()-0.5) * radius

    return x, y

def polynomialPotential(value,firstDegreeFactor=0.0,secondDegreeFactor=0.0,thirdDegreeFactor=0.0):
    retval = 0.0

    if type(firstDegreeFactor) == float and firstDegreeFactor != 0.0:
        retval += firstDegreeFactor * value

    if type(secondDegreeFactor) == float and secondDegreeFactor != 0.0:
        retval += secondDegreeFactor * value ** 2

    if type(thirdDegreeFactor) == float and thirdDegreeFactor != 0.0:
        retval += thirdDegreeFactor * value ** 3

    return retval





def testGraph(testType):

    nw = NodeWrangler()
    nw.genNodes(40,10)

    origionalNodes = copy.deepcopy(nw.nodes)

    nw.genAllPhysicalLinks(1.0)
    print(nw)
    #nw.CalculateForceVectors(1.0)
    #nw.ApplyVectors()

    nw.CalculateForceVectors(10.0)

    #nw.iterativeCalculateForceVectors(1,0.01,2.00)
    #nw.iterativeCalculateForceVectors(1,0.5,0.01)

    #print(nw)

    nw.dispStreamLineGraph()




print("Some trig testing:")
n = 100

for i in range(n):
    
    


