import math
import random
import copy

import NEWLIFE_IDtrack as IDtrack

import matplotlib.pyplot as pplot
import matplotlib.backends.backend_qt5agg
import numpy

import SomeExternalLibs.graphics as gfx


#from sortedcontainers import SortedDict



class NodeWrangler(object):
    def __init__(self):

        self.nodes = IDtrack.IDtracker()
        self.links = IDtrack.IDtracker()

        self._originalNodes = None

        self.fig = None

    def __str__(self):
        return "\n".join(str(self.nodes[n]) for n in self.nodes)

    def genNodes(self,n:int,radius:float):
        for i in range(n):
            rn = genRandInsideUnitCircle(radius)
            self.nodes.StoreVal(PhysicalNode(rn, 1.0, None))
        self._originalNodes = copy.deepcopy(self.nodes)

    def genAllChargedLinks(self, strength:float):

        #This guy can be optimized eventually by making a todo list that is iterated through and have values poped as the list is iterated through.


        completedkeys = []

        for key in self.nodes:

            completedkeys.append(key)
            connectToList = {k: self.nodes[k] for k in set(self.nodes.keys()) - set(completedkeys)}

            for connection in connectToList:

                pl = physicalNodeLink(self.nodes[key],self.nodes[connection])

                linkkey = self.links.StoreVal(pl)

                self.nodes[key].links.append(linkkey)
                self.nodes[connection].links.append(linkkey)

    def removeAllLinksOfInstance(self, instanceType: type):

        #literally black magic. Do not touch!

        l = list(self.links.keys())
        iternum = 0
        indexPicker = 0
        selection=0

        while indexPicker <= len(self.links)-1-selection:
            selection = iternum-indexPicker
            if isinstance(self.links[iternum], instanceType):

                #print(self.links[selection].n1.links)
                self.links[selection].n1.links.pop(self.links[selection].n1.links.index(selection))
                self.links[selection].n2.links.pop(self.links[selection].n2.links.index(selection))

                self.links.delVal(selection)

            else:
                indexPicker+=1
            iternum+=1

    def removeAllLinksOfType(self, instanceType: type):

        #Also Black Magic

        l = list(self.links.keys())
        iternum = 0
        indexPicker = 0
        selection=0

        while indexPicker <= len(self.links)-1-selection:

            selection = iternum - indexPicker

            if type(self.links[iternum]) == instanceType:

                #print(self.links[selection].n1.links)
                self.links[selection].n1.links.pop(self.links[selection].n1.links.index(selection))
                self.links[selection].n2.links.pop(self.links[selection].n2.links.index(selection))

                self.links.delVal(selection)

            else:
                indexPicker+=1
            iternum+=1

    def get_unlinked_nodes(self):
        return [n for n in self.nodes if len(self.nodes[n].links) == 0]




    def genClosestNSpringLinks(self,N):
        '''
        try:

            from sortedcollections import sortedDict

            for n in self.nodes:

                print("n={}/{}\n".format(n,len(self.nodes)-1))
                lenDict = {}
                for l in self.nodes[n].links:
                    lenDict[l]= self.links[l].length()

                #k: self.links[self.nodes[n].links[k]].length()for k in
                sortedShortestKeys = SortedDict(lenDict).keys()[0:N]

            print("linksd len is {}".format(len(linksd)))

            #lengths="\n\t{}".join("" for j in linksd).format([lenDict[j] for j in linksd])


            print("there are {} links their lengths are :".format(len(sortedShortestKeys)))#,lengths))
            for i in sortedShortestKeys:
                print(i,", ",lenDict[i])



        except:
        '''
        for n in self.nodes:

            print("n={}/{}total nodes. There are {} links".format(n, len(self.nodes)-1,len(self.nodes[n].links)))

            lenDict = {k:self.links[k].lengthSquared() for k in self.nodes[n].links}

            #vallist = [val for val in lenDict.values()]

            LowestVal = None
            LowestKeys = []
            
            for key in lenDict:
                if len(LowestKeys) == N:
                    break

                if LowestVal is None:
                    LowestVal = lenDict[key]

                elif key not in LowestKeys:
                    if lenDict[key] < LowestVal:
                        LowestVal = lenDict[key]
                        LowestKeys.append(key)


            print("{} keys were selected at node {}".format(len(LowestKeys),n))
            for k in LowestKeys:

                print("The selected key was {} with a squared length of {}\n".format(k,lenDict[k]))

                self.links[k] = SpringLink(self.links[k].n1, self.links[k].n2, desiredSpringLength=1.0, springStrength=1.0)

                #self.nodes[n].links.append()

    @classmethod
    def angFromSinCos(self, sinang:float, cosang:float):

        #sinang,cosang = math.sin(sincomp),math.cos(coscomp)

        # Quadrant breakup
        #Q1
        if sinang > 0.0 and cosang >= 0.0:
            #print("Q1")
            return math.cos(cosang)
        # Q2
        elif sinang > 0.0 and cosang <= 0.0:
            #print("Q2")
            return math.pi - math.cos(cosang)
        # Q3
        elif sinang < 0.0 and cosang <= 0.0:
            #print("Q3")
            return math.pi + math.cos(cosang)
        # Q4
        elif sinang < 0.0 and cosang >= 0.0:
            #print("Q4")
            return 2 * math.pi - math.cos(cosang)

        else:
            print("*Error with values of sin and cos components: {} {}*".format(sinang,cosang))
            raise ValueError

    def resetMoveVectors(self):
        for n in self.nodes:
            #print("Node: {}".format(self.nodes[n].id))
            #print("movVec: {}".format(type(self.nodes[n].moveVec)))
            self.nodes[n].moveVec = Pos2D((0, 0))

    def iterativeCalculateForceVectors(self,iterations:int,startStr:float,endStr:float):

        for x in range(iterations):
            m=(endStr-startStr)/iterations
            b=startStr

            self.CalculateForceVectors(m*x+b)
            #print("m,x,b = {},{},{}".format(m,x,b))
            self.ApplyVectors()


    def CalculateForceVectors(self, strength):
        print("Iterating through vectors for force calculation")

        for l in self.links:

            #(l.n1 * l.n2)/l.lengthSquared()

            hypSq = self.links[l].lengthSquared()
            hyp = math.sqrt(hypSq)

            sincomp = self.links[l].deltaY / hyp
            coscomp = self.links[l].deltaX / hyp
            #print(sincomp, coscomp)

            ang = NodeWrangler.angFromSinCos(sincomp, coscomp)

            moveMag = 0.0

            # need to recalculate hypsquared here so optional parameter to ease calculation
            if isinstance(self.links[l],physicalNodeLink):
                moveMag += self.links[l].chargedParticleForce(hypSq)

            if isinstance(self.links[l],SpringLink):
                #print("spring strength applied")
                moveMag += self.links[l].getSpringStrength(hyp)



            #springForce = self.links[l]




            #print("cosang: {} sinang: {} chargeMul: {} ang: {}".format(coscomp*180/math.pi,sincomp*180/math.pi,chargemul,ang*180/math.pi))


            #print("N1:{},N2:{} ".format(self.links[l].n1.moveVec,self.links[l].n2.moveVec))

            self.links[l].n1.moveVec.x += -moveMag * strength / self.links[l].n1.weight * math.cos(ang)
            self.links[l].n1.moveVec.y += -moveMag * strength / self.links[l].n1.weight * math.sin(ang)

            self.links[l].n2.moveVec.x += moveMag * strength / self.links[l].n2.weight * math.cos(ang)
            self.links[l].n2.moveVec.y += moveMag * strength / self.links[l].n2.weight * math.sin(ang)


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
        data = self.genQuiverGraphData()
        pplot.scatter(data[0], data[1])
        pplot.show()

    def dispQuiver(self):
        data = self.genQuiverGraphData()
        pplot.quiver(data[0], data[1], data[2], data[3], label = "This label")
        pplot.show()

    def assignMoveVecsToDisplacementFromOrig(self):

        self.ApplyVectors()

        for n in self.nodes:
            print("n={}".format(n))

            self.nodes[n].moveVec = Pos2D((self.nodes[n].x-self._originalNodes[n].x, self.nodes[n].y-self._originalNodes[n].y))


    def genQuiverGraphData(self):
        xdata = [self.nodes[n].x for n in self.nodes]
        ydata = [self.nodes[n].y for n in self.nodes]

        xmoveData = [self.nodes[n].moveVec.x for n in self.nodes]
        ymoveData = [self.nodes[n].moveVec.y for n in self.nodes]

        data = numpy.array([xdata, ydata,xmoveData,ymoveData])

        print("some graph data:\n\t{},{} pieces of data in xdata,ydata\n\t{},{} pieces of data in xmovedata,ymovedata".format(len(xdata),len(ydata),len(xmoveData),len(ymoveData)))

        #p = pplot.quiver(data[0], data[1], data[2], data[3])
        return data

class Plotter(object):
    def __init__(self, nwo: NodeWrangler):
        self.can = gfx.GraphWin("WorldGen", 1000, 800)

        # can.mainloop()

        nodestoplot = {k: [v.pos,v.weight] for k, v in nwo.nodes.items()}
        self.plotPoints(nodestoplot, 20)

        #print(nwo.nodes.keys())

        linestoplot = {k: [v.n1, v.n2] for k, v in nwo.links.items()}
        self.plotLines(linestoplot, 20)

        self.can.getMouse()

    def plotPoints(self,point_id_dict: dict, scale= 1.0):

        for n in point_id_dict:
            c = gfx.Circle(center=gfx.Point(scale*point_id_dict[n][0][0]+self.can.width//2, scale*point_id_dict[n][0][1]+self.can.height//2), radius=scale/5*point_id_dict[n][1])
            c.draw(self.can)

    def plotLines(self, line_id_dict: dict, scale=1.0):

        for l in line_id_dict:
            #print(point_id_dict[0][0])
            line = gfx.Line(gfx.Point(scale*line_id_dict[l][0][0]+self.can.width//2, scale*line_id_dict[l][0][1]+self.can.height//2),gfx.Point(scale*line_id_dict[l][1][0]+self.can.width//2, scale*line_id_dict[l][1][1]+self.can.height//2))

            line.draw(self.can)


#Pos2D class Collection
class Pos2D(object):

    def __init__(self,pos):
        #self.pos = None
        self.__pos = None

        #This is a property and sets __pos
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

#NodeLink class Collection
class nodeLink(object):

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
    def __init__(self,n1, n2):

        super().__init__(n1, n2)

    def chargedParticleForce(self, lenSquared:float = None):

        if lenSquared is not None:
            return (self.n1.weight * self.n2.weight) / lenSquared
        else:
            return (self.n1.weight * self.n2.weight) / self.lengthSquared()

class SpringLink(physicalNodeLink):

    def __init__(self,n1:PhysicalNode,n2:PhysicalNode,desiredSpringLength:float,springStrength:float):
        super().__init__(n1,n2)
        self.desiredLength = desiredSpringLength
        self.springStrength = springStrength

    def springDisplacement(self, length:float = None):
        if length is None:
            return self.length() - self.desiredLength
        else:
            return self.desiredLength - length

    def getSpringStrength(self,length:float=None):
        return self.springDisplacement(length) * self.springStrength





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

def testGraph():
    print(pplot.get_backend())

    nw = NodeWrangler()
    nw.genNodes(50, 10)
    nw.genAllChargedLinks(1.0)

    #print(nw)

    #nw.CalculateForceVectors(1.0)
    #nw.ApplyVectors()


    #nw.CalculateForceVectors(10.0)





    nw.genClosestNSpringLinks(1)


    nw.iterativeCalculateForceVectors(10, 0.0, 0.5)
    nw.iterativeCalculateForceVectors(10, 0.5, 0.0)


    nw.removeAllLinksOfType(physicalNodeLink)

    print("there are {} links".format(len(nw.links)))

    print(nw.get_unlinked_nodes())

    #nw.assignMoveVecsToDisplacementFromOrig()
    #nw.dispQuiver()
    #nw.resetMoveVectors()

    Plotter(nw)

    #nw.dispQuiver()
    #nw.fig.show()



testGraph()
