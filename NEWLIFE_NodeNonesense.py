import math
import random
import copy

import NEWLIFE_IDtrack as IDtrack

import matplotlib.pyplot as pplot
import matplotlib.backends.backend_qt5agg
import numpy

import SomeExternalLibs.graphics as gfx

class NodeWrangler(object):
    def __init__(self):

        self.nodes = IDtrack.IDtracker()
        self.links = IDtrack.IDtracker()

        self._originalNodes = None

        self.fig = None

    def __str__(self):
        return "\n".join(str(self.nodes[n]) for n in self.nodes)

    def gen_nodes(self, n:int, radius:float):
        for i in range(n):
            rn = gen_rand_pos_in_unit_circle(radius)
            self.nodes.StoreVal(PhysicalNode(rn, 1.0, None))
        self._originalNodes = copy.deepcopy(self.nodes)


    #
    def gen_physical_links_between_all(self, strength:float):

        #todo This guy can be optimized eventually by making a list that is iterated through and have values poped as the list is iterated through.


        completedkeys = []

        for key in self.nodes:

            completedkeys.append(key)
            connectToList = {k: self.nodes[k] for k in set(self.nodes.keys()) - set(completedkeys)}

            for connection in connectToList:

                pl = PhysicalNodeLink(self.nodes[key], self.nodes[connection])

                linkkey = self.links.StoreVal(pl)

                self.nodes[key].links.append(linkkey)
                self.nodes[connection].links.append(linkkey)

    def remove_links_of_instance(self, instanceType: type):

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

    def remove_links_of_type(self, instanceType: type):

        #Also Black Magic

        l = list(self.links.keys())
        iternum = 0
        indexPicker = 0
        selection=0

        while indexPicker <= len(self.links)-1:

            selection = iternum - indexPicker

            if type(self.links[iternum]) == instanceType:

                #print(self.links[selection].n1.links)
                self.links[selection].n1.links.pop(self.links[selection].n1.links.index(selection))
                self.links[selection].n2.links.pop(self.links[selection].n2.links.index(selection))

                self.links.delVal(selection)

            else:
                indexPicker+=1
            iternum+=1

    #Returns a list of all nodes that are not linked
    def get_unlinked_nodes(self):
        return [n for n in self.nodes if len(self.nodes[n].links) == 0]

    # Returns a list of the number of links of each node
    def get_node_links(self):
        return [len(self.nodes[n].links) for n in self.nodes]

    def gen_closest_n_spring_links(self, N):

        for n in self.nodes:

            print("n={}/{}total nodes. There are {} links".format(n, len(self.nodes)-1,len(self.nodes[n].links)))


            #creates a dictionary of connections of the current node as LinkKey: LinkLength
            lenDict = {k:self.links[k].lengthSquared() for k in self.nodes[n].links}

            #LowestVal is the current lowest length found in the iteration
            lowest_key = None

            #LowestKeys is the list of keys that have been determined to be the lowest
            LowestKeys = []

            #Loop until the correct number of lengths have been found
            while len(LowestKeys) != N:

                #Iterate through the keys of the node links
                for key in lenDict:

                    #Assigns the lowest val upon the first iteration
                    if lowest_key is None:
                        lowest_key = key

                    #Check that
                    elif lenDict[key] < lenDict[lowest_key]:

                        #Check that current key hasn't been used as a previous lowest
                        if key not in LowestKeys:

                            #assign a new lowest key
                            lowest_key = key

                #After the entire list has been iterated over, add the lowest value to the list of lowest keys
                LowestKeys.append(lowest_key)

            print("{} keys were selected at node {}".format(len(LowestKeys),n))
            for k in LowestKeys:

                print("The selected key was {}, connected between {} and {} with a squared length of {}\n".format(k,self.links[k].n1,self.links[k].n2,lenDict[k]))
                #self.links.delVal(k)

                #temp = self.links[k]

                self.links[k] = SpringLink(self.links[k].n1, self.links[k].n2, desiredSpringLength=1.0, springStrength=1.0)

                #self.nodes[n].links.append()

    @classmethod
    def ang_from_sin_cos_ratios(self, sinang:float, cosang:float):

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

    def reset_move_vectors(self):
        for n in self.nodes:
            self.nodes[n].moveVec = Pos2D((0.0, 0.0))


    def iterative_calculate_force_vectors(self, iterations:int, startStr:float, endStr:float):
        '''
        Iteratively applies forces to nodes based on connections

        :param iterations:
        :param startStr: strength of force application at the start of the iteration
        :param endStr: strength of force application at the end of the iteration
        :return:
        '''

        for x in range(iterations):
            m=(endStr-startStr)/iterations
            b=startStr

            self.calc_force_vectors(m * x + b, m * x + b)
            self.apply_vectors()


    def calc_force_vectors(self, charge_strength, spring_strength):
        print("Iterating through vectors for force calculation")

        for l in self.links:

            #(l.n1 * l.n2)/l.lengthSquared()

            hypSq = self.links[l].lengthSquared()
            hyp = math.sqrt(hypSq)

            sincomp = self.links[l].deltaY / hyp
            coscomp = self.links[l].deltaX / hyp
            #print(sincomp, coscomp)

            ang = NodeWrangler.ang_from_sin_cos_ratios(sincomp, coscomp)

            moveMag = 0.0

            # need to recalculate hypsquared here so optional parameter to ease calculation
            if isinstance(self.links[l], PhysicalNodeLink):
                moveMag += charge_strength * self.links[l].charged_particle_force(hypSq)

            if isinstance(self.links[l],SpringLink):
                #print("spring strength applied")
                moveMag += spring_strength * self.links[l].get_spring_strength(hyp)

            #print("cosang: {} sinang: {} chargeMul: {} ang: {}".format(coscomp*180/math.pi,sincomp*180/math.pi,chargemul,ang*180/math.pi))
            #print("N1:{},N2:{} ".format(self.links[l].n1.moveVec,self.links[l].n2.moveVec))

            self.links[l].n1.moveVec.x += -moveMag / self.links[l].n1.weight * math.cos(ang)
            self.links[l].n1.moveVec.y += -moveMag / self.links[l].n1.weight * math.sin(ang)

            self.links[l].n2.moveVec.x += moveMag / self.links[l].n2.weight * math.cos(ang)
            self.links[l].n2.moveVec.y += moveMag / self.links[l].n2.weight * math.sin(ang)

    def apply_vectors(self):

        #Apply the moveVectors
        for n in self.nodes:

            if self.nodes[n].rooted == False:

                #print("type: {}, node x: {}, y: {}".format(type(self.nodes[n]),self.nodes[n].x, self.nodes[n].y)) #+= self.nodes[n].moveVec.x
                #print("type: {}, movevec x: {}, y: {}".format(type(self.nodes[n].moveVec),self.nodes[n].moveVec.x, self.nodes[n].moveVec.y))#+= self.nodes[n].moveVec.y
                self.nodes[n].x += self.nodes[n].moveVec.x
                self.nodes[n].y += self.nodes[n].moveVec.y

        self.reset_move_vectors()

    def disp_plot(self):
        data = self.gen_quiver_graph_data()
        pplot.scatter(data[0], data[1])
        pplot.show()

    def disp_quiver(self):
        data = self.gen_quiver_graph_data()
        pplot.quiver(data[0], data[1], data[2], data[3], label = "This label")
        pplot.show()

    def assign_moveVec_to_displacement_from_origional_position(self):

        self.apply_vectors()
        for n in self.nodes:
            print("n={}".format(n))
            self.nodes[n].moveVec = Pos2D((self.nodes[n].x-self._originalNodes[n].x, self.nodes[n].y-self._originalNodes[n].y))

    def gen_quiver_graph_data(self):
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
        self.plot_points(nodestoplot, 20)

        #print(nwo.nodes.keys())

        linestoplot = {k: [v.n1, v.n2] for k, v in nwo.links.items()}
        self.plot_lines(linestoplot, 20)

        self.can.getMouse()

    def plot_points(self, point_id_dict: dict, scale= 1.0):

        for n in point_id_dict:
            c = gfx.Circle(center=gfx.Point(scale*point_id_dict[n][0][0]+self.can.width//2, scale*point_id_dict[n][0][1]+self.can.height//2), radius=scale/5*point_id_dict[n][1])
            c.draw(self.can)

    def plot_lines(self, line_id_dict: dict, scale=1.0):

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
        if type(id) is None:
            return str(super())
        else:
            return "id: {}, x: {}, y: {}".format(self.id, self.x, self.y)

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

class PhysicalNodeLink(nodeLink):
    def __init__(self,n1, n2):

        super().__init__(n1, n2)

    def charged_particle_force(self, lenSquared:float = None):

        if lenSquared is not None:
            return (self.n1.weight * self.n2.weight) / lenSquared
        else:
            return (self.n1.weight * self.n2.weight) / self.lengthSquared()

class SpringLink(PhysicalNodeLink):

    def __init__(self,n1:PhysicalNode,n2:PhysicalNode,desiredSpringLength:float,springStrength:float):
        super().__init__(n1,n2)
        self.desiredLength = desiredSpringLength
        self.springStrength = springStrength

    def spring_displacement(self, length:float = None):
        if length is None:
            return self.length() - self.desiredLength
        else:
            return self.desiredLength - length

    def get_spring_strength(self, length:float=None):
        return self.spring_displacement(length) * self.springStrength


def gen_rand_pos_in_unit_circle(radius=1.0):
    x = (random.random() - 0.5) * radius
    y = (random.random() - 0.5) * radius

    while x ** 2 + y ** 2 > radius ** 2:
        x = (random.random()-0.5) * radius
        y = (random.random()-0.5) * radius

    return x, y

def polynomialPotential(value,bias:float = 0.0,firstDegreeFactor:float=0.0,secondDegreeFactor:float=0.0,thirdDegreeFactor:float=0.0):

    retval = bias

    if type(firstDegreeFactor) == float and firstDegreeFactor != 0.0:
        retval += firstDegreeFactor * value

    if type(secondDegreeFactor) == float and secondDegreeFactor != 0.0:
        retval += secondDegreeFactor * value ** 2

    if type(thirdDegreeFactor) == float and thirdDegreeFactor != 0.0:
        retval += thirdDegreeFactor * value ** 3

    return retval

def test_graph():

    nw = NodeWrangler()
    nw.gen_nodes(5, 10)
    nw.gen_physical_links_between_all(1.0)

    nw.gen_closest_n_spring_links(1)

    #nw.iterative_calculate_force_vectors(10, 0.0, 0.5)
    #nw.iterative_calculate_force_vectors(10, 0.5, 0.0)

    nw.remove_links_of_type(PhysicalNodeLink)

    print("there are {} links".format(len(nw.links)))

    print(nw.get_unlinked_nodes())
    print(nw.get_node_links())

    Plotter(nw)

def plot_quiver():
    print(pplot.get_backend())

    nw = NodeWrangler()
    nw.gen_nodes(5, 10)
    nw.gen_physical_links_between_all(1.0)

    nw.iterative_calculate_force_vectors(10, 0.0, 0.5)
    nw.iterative_calculate_force_vectors(10, 0.5, 0.0)


    nw.assign_moveVec_to_displacement_from_origional_position()
    nw.disp_quiver()
    nw.reset_move_vectors()

if __name__ == "__main__":
    test_graph()



"""
hi reddit.

The bit I am having trouble with is on line 111: gen_closest_n_spring_links()

The code based on the debugger seems to be finding the correct link.

But it is not creating the link between the nodes it selects

Furthermore, the first node never generates a link, and the second only sometimes does

The later links also generate more than a single link.

I am overlooking something HUGE but can not put my finger on it. This has taken me the better part of the afternoon and I would really appreciate another pair of eyes.

General code review is also welcomed.


Here is the output I am getting:
    
    n=0/4total nodes. There are 4 links
    1 keys were selected at node 0
    The selected key was 3, connected between id: 0, x: 3.4777810915689225, y: -2.9876675771294403 and id: 4, x: 2.7149802188741834, y: -2.189168706378423 with a squared length of 1.2194656179745056
    
    n=1/4total nodes. There are 4 links
    1 keys were selected at node 1
    The selected key was 5, connected between id: 1, x: -2.815756580592086, y: 3.398505504019141 and id: 3, x: -0.7898316784048398, y: 3.439405419356868 with a squared length of 4.106044512377037
    
    n=2/4total nodes. There are 4 links
    1 keys were selected at node 2
    The selected key was 4, connected between id: 1, x: -2.815756580592086, y: 3.398505504019141 and id: 2, x: -3.9549320995259674, y: -1.7672549140054683 with a squared length of 27.982801559368063
    
    n=3/4total nodes. There are 4 links
    1 keys were selected at node 3
    The selected key was 5, connected between id: 1, x: -2.815756580592086, y: 3.398505504019141 and id: 3, x: -0.7898316784048398, y: 3.439405419356868 with a squared length of 4.106044512377037
    
    n=4/4total nodes. There are 4 links
    1 keys were selected at node 4
    The selected key was 3, connected between id: 0, x: 3.4777810915689225, y: -2.9876675771294403 and id: 4, x: 2.7149802188741834, y: -2.189168706378423 with a squared length of 1.2194656179745056
    
    there are 3 links
    [0, 1]
    [0, 0, 2, 2, 2]

"""