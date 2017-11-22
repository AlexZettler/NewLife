import math
import random
import copy
import abc

import IDtrack as IDtrack

import matplotlib.pyplot as pplot
import numpy

import SomeExternalLibs.graphics as gfx


class NodeWrangler(object):
    def __init__(self):

        self.nodes = IDtrack.id_tracker()
        self.links = IDtrack.id_tracker()

        self._originalNodes = None

        #self.fig = None

    def __str__(self):
        return "\n".join(str(self.nodes[n]) for n in self.nodes)

    def getStats(self):

        print("Getting you some stats of the graph...\n\t" + "\n\t".join(
            ("there are {} links".format(len(self.links)),
             "unlinked nodes: {}".format(self.get_unlinked_nodes()),
             "number of links per node: {}".format(self.get_node_links()),
             "number of islands: {}".format(self.get_node_islands()))))


    def gen_nodes_with_radius(self, n:int, radius:float):
        for i in range(n):
            rn = Pos2D.gen_rand_pos_in_unit_circle(radius)
            self.nodes.store_val(PhysicalNode(rn.x, rn.y, 1.0))
        self._originalNodes = copy.deepcopy(self.nodes)

    def gen_nodes_with_params(self, n=None, radius=None, density=None):

        if n is None:

            if radius is None or density is None:
                raise ValueError
            else:
                self.gen_nodes_with_radius(density / (math.pi * radius ** 2), radius)

        elif radius is None:

            if n is None or density is None:
                raise ValueError
            else:

                (math.pi * radius ** 2)

                self.gen_nodes_with_radius(n, )

        elif density is None:

            if n is None or radius is None:
                raise ValueError
            else:
                self.gen_nodes_with_radius(n, radius)

        else:
            raise ValueError

    @classmethod
    def ge_links_between_all(self, tracker, link_type, **kwargs):

        #todo This guy can be optimized eventually by making a list that is iterated through and have values poped as the list is iterated through.

        completedkeys = []

        if isinstance(link_type(Pos2D(0, 0), Pos2D(0, 0), **kwargs), NodeLink):

            for key in tracker.nodes:

                completedkeys.append(key)

                connectToList = {k: tracker.nodes[k] for k in set(tracker.nodes) - set(completedkeys)}

                for connection_key in connectToList:
                    pl = link_type(tracker.nodes[key], tracker.nodes[connection_key], **kwargs)

                    tracker.links.store_val(pl)

                    # print("pl={}".format(pl.get_id()))
                    tracker.nodes[key].links.append(pl.get_id())
                    tracker.nodes[connection_key].links.append(pl.get_id())

        else:
            raise ValueError

    def remove_links_of_type(self, instance_type: type):
        keysToRemove=[]

        for l in self.links:
            if type(self.links[l]) is instance_type:
                keysToRemove.append(l)
                #print(self.links[selection].n1.links)

        for l in keysToRemove:
            self.remove_link(l)

    def remove_link(self,link_id):
        self.links[link_id].n1.links.pop(self.links[link_id].n1.links.index(link_id))
        self.links[link_id].n2.links.pop(self.links[link_id].n2.links.index(link_id))
        self.links.del_val(link_id)

    #Returns a list of all nodes that are not linked
    def get_unlinked_nodes(self):
        return [n for n in self.nodes if len(self.nodes[n].links) == 0]

    # Returns a list of the number of links of each node
    def get_node_links(self):
        return {n:len(self.nodes[n].links) for n in self.nodes}

    def get_node_islands(self):

        '''
        :return: A dictionary containing a list of nodes for each island key
        '''

        open_set = []
        closed_set = []

        island_counter = 0

        # a list of length[number of nodes] that maps nodes to an island group
        node_islands = {}

        for n in self.nodes:

            if n not in closed_set:

                open_set.insert(0, n)
                node_islands[island_counter] = []

                while len(open_set) > 0:

                    open_set_selection = open_set.pop()
                    # print("ossel={}".format(open_set_selection))

                    for l in self.nodes[int(open_set_selection)].links:

                        # print("nodeID = {}, linkID = {}".format(open_set_selection, l))

                        c = self.links[l].get_connected(open_set_selection).get_id()

                        if c not in closed_set:

                            open_set.insert(0, c)
                            node_islands[island_counter].append(int(c))
                        closed_set.append(c)

                island_counter += 1

        #node_island_links = IDtrack.IDtracker()

        return node_islands

    def connect_node_islands(self):

        while True:

            ni = self.get_node_islands()

            if len(ni) <= 1:
                break


            island_wrangler = NodeWrangler()

            # create nodegroups with collection of nodes from each island
            for i in ni:
                island_wrangler.nodes.store_val(NodeGroup(self, [n for n in ni[i]], None))


            # create a lose link to find shortest connections
            island_wrangler.ge_links_between_all(island_wrangler, NodeLink)

            #connect all closest node island groups
            for island_node in island_wrangler.nodes:
                island_wrangler.connect_links(
                    island_wrangler.gen_closest_n_links(1, island_wrangler.nodes, island_wrangler.links), SpringLink,
                    springStrength=0.0, desiredLength=0.0)

            #
            island_wrangler.remove_links_of_type(NodeLink)

            # Plotter(island_wrangler, 15)

            # this handles indexing of changed keys
            node_group_mapping = {k: island_wrangler.nodes[k].get_id() for k in island_wrangler.nodes}

            for ng in island_wrangler.nodes:

                for node_group_link in island_wrangler.nodes[ng].links:

                    first_group = node_group_mapping[ng]
                    second_group = node_group_mapping[island_wrangler.links[node_group_link].get_connected(ng).get_id()]

                    if first_group == second_group:
                        print("\n")
                        print("node group={}, node group link={}".format(ng, node_group_link))

                        print("{} indexed for list of nodes: {}".format(island_wrangler.nodes[ng].get_id(),
                                                                        island_wrangler.nodes))

                        print("{} indexed for list of links: {}".format(island_wrangler.links[node_group_link],
                                                                        island_wrangler.nodes[ng].links))

                    if first_group != second_group:
                        island_wrangler.nodes[ng].true_closet_connect(self, island_wrangler.nodes[first_group],
                                                                      island_wrangler.nodes[second_group])
                        node_group_mapping[second_group] = first_group


            island_wrangler.remove_links_of_type(PhysicalNodeLink)

            Plotter(island_wrangler, 15)

                #se = NodeGroup(self,node_islands[k])

    #@classmethod
    def gen_closest_n_links(self, N, node_tracker, link_tracker):

        retlinks = []

        for n in self.nodes:

            print("n={}/{}total nodes. There are {} links".format(n, len(self.nodes)-1,len(self.nodes[n].links)))

            #creates a dictionary of connections of the current node as LinkKey: LinkLength
            lenDict = {k:self.links[k].lengthSquared() for k in self.nodes[n].links}

            lowest_key=None

            #LowestKeys is the list of keys that have been determined to be the lowest
            LowestKeys = []

            #Loop until the correct number of lengths have been found
            while len(LowestKeys) != N or lowest_key is None:

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

            #print("{} keys were selected at node {}".format(len(LowestKeys),n))
            for k in LowestKeys:
                #print("The selected key was {}, connected between {} and {} with a squared length of {}\n".format(k,self.links[k].n1,self.links[k].n2,lenDict[k]))

                #if type(self.links[k]) == PhysicalNodeLink:
                retlinks.append(k)

        return retlinks

    def connect_links(self, links:list, link_type:type, **kwargs):

        if isinstance(link_type(Pos2D(0, 0), Pos2D(0, 0), **kwargs), NodeLink):

            for l in links:

                self.links[l] = link_type(self.links[l].n1, self.links[l].n2, **kwargs)

        else:

            raise TypeError

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
            self.nodes[n].moveVec = Pos2D(0.0, 0.0)

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

            if hypSq == 0.0:
                print("n1={}, n2={}".format(self.links[l].n1, self.links[l].n2))


            hyp = math.sqrt(hypSq)

            #try:

            sincomp = self.links[l].deltaY / hyp
            coscomp = self.links[l].deltaX / hyp

            # except ZeroDivisionError:

            #print("something happened")

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
    def __init__(self, nwo: NodeWrangler, scale):
        self.can = gfx.GraphWin("WorldGen", 1000, 800)

        # can.mainloop()

        if isinstance(nwo.nodes[0], PhysicalNode):
            nodestoplot = {k: [v.pos, v.weight] for k, v in nwo.nodes.items()}
        elif isinstance(nwo.nodes[0], NodeGroup):
            nodestoplot = {k: [v.pos, 1.0] for k, v in nwo.nodes.items()}
        else:
            raise TypeError

        self.plot_points(nodestoplot, scale)

        #print(nwo.nodes.keys())

        linestoplot = {k: [v.n1, v.n2] for k, v in nwo.links.items()}
        self.plot_lines(linestoplot, scale)

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
    def __init__(self, x, y):
        self.x, self.y = float(x), float(y)

    def __str__(self):
        return "Pos2D-x:{}, y:{}".format(self.x, self.y)
        pass

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        else:
            raise ValueError

    def __setitem__(self, key, value):
        if type(key) is int:
            if key == 0:
                self.x = value
            elif key == 1:
                self.y = value
        else:
            raise ValueError

    def __iadd__(self, other):
        if isinstance(other, Pos2D):
            self.x += other.x
            self.y += other.y
            return self

    def __isub__(self, other):
        if isinstance(other, Pos2D):
            self.x -= other.x
            self.y -= other.y
            return self
        else:
            raise TypeError

    def __itruediv__(self, other):
        if isinstance(other, Pos2D):
            self.x /= other.x
            self.y /= other.y
        elif type(other) is float:
            self.x /= other
            self.y /= other
        else:
            raise TypeError
        return self

    @property
    def pos(self):
        return (self.x, self.y)

    @pos.setter
    def pos(self, val):
        if len(val) == 2 and type(val) in (list, tuple):
            self.x, self.y = val[0], val[1]

        else:
            raise ValueError

    def __len__(self):
        return 2

    @classmethod
    def gen_rand_pos_in_unit_circle(cls, radius=1.0):

        x = (random.random() - 0.5) * radius
        y = (random.random() - 0.5) * radius

        while x ** 2 + y ** 2 > radius ** 2:
            x = (random.random() - 0.5) * radius
            y = (random.random() - 0.5) * radius

        return Pos2D(x, y)


class Node(Pos2D, IDtrack.tracked_object):
    def __init__(self, x, y):

        super().__init__(x,y)

        self.links = []

    def track(self, id):
        IDtrack.tracked_object.__init__(self, id)

    def __str__(self):
        if id is None:
            return super().__str__()
        else:
            return "id: {}, x: {}, y: {}".format(self.get_id(), self.x, self.y)

    def __int__(self):
        return self.get_id()

class PhysicalNode(Node):
    def __init__(self, x, y, weight):
        Node.__init__(self, x, y)

        self.weight = weight
        self.moveVec = Pos2D(0.0, 0.0)
        self.rooted = False

    def __str(self):
        return "{}, moveVec: {}".format(str(super), self.moveVec)


class NodeGroup(Node):
    # todo Build this as a group of nodes acting as a single point in the center of the nodes. The probelm will need to be solved by finding nodes closest to center of another NodeGroup and connecting it to the
    # todo Debug
    def __init__(self, wranglerParent: NodeWrangler, Nodes: list, id):

        self.Nodes = Nodes
        center = Pos2D(0.0, 0.0)

        for k in self.Nodes:
            center += wranglerParent.nodes[k]
        # print(len(Nodes))
        center /= float(len(Nodes))

        Node.__init__(self, center.x, center.y)

    def true_closet_connect(self, nodeWranglerParent: NodeWrangler, first_group, second_group):

        shortest_link = None

        NL = None

        print(first_group.Nodes)
        print(second_group.Nodes)


        #iterate through all possible node connections
        for n1 in first_group.Nodes:

            for n2 in second_group.Nodes:

                # print("n1={},n2={}".format(n1, n2))

                NL = NodeLink(nodeWranglerParent.nodes[n2], nodeWranglerParent.nodes[n1])

                if shortest_link is None:
                    shortest_link = NL

                elif NL.lengthSquared() < shortest_link.lengthSquared():

                    shortest_link = NL

        if NL is not None:

            NL = SpringLink(NL.n1, NL.n2, desiredLength=1.0, springStrength= 1.0)

            nodeWranglerParent.links.store_val(NL)
            nodeWranglerParent.nodes[NL.n1.get_id()].links.append(NL.get_id())
            nodeWranglerParent.nodes[NL.n2.get_id()].links.append(NL.get_id())

            first_group.Nodes += second_group.Nodes

        else:
            raise ValueError


#NodeLink class Collection
class NodeLink(IDtrack.tracked_object):

    def __init__(self,n1,n2,**kwargs):
        self.n1, self.n2 = n1, n2

        if n1 == n2:
            raise ValueError

        IDtrack.tracked_object.__init__(self, None)

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

    def get_connected(self, n):
        intn1 = int(self.n1)
        intn2 = int(self.n2)
        intn = int(n)

        # print("Trying to get connection from {}, between {} and {}".format(intn,intn1, intn2))

        # print(intn == intn1)

        if intn == intn1:
            retid = self.n2
        elif intn == intn2:
            retid = self.n1
        else:

            raise KeyError

        #print("Printing the connected link of node {}\n\tThe connected node is: {}\n\tThe link id is: {}".format(self.n1, self.n2, retid))
        return retid


class PhysicalNodeLink(NodeLink):
    def __init__(self,n1, n2, **kwargs):

        super().__init__(n1, n2, **kwargs)

    def charged_particle_force(self, lenSquared:float = None):

        if lenSquared is not None:
            return (self.n1.weight * self.n2.weight) / lenSquared
        else:
            return (self.n1.weight * self.n2.weight) / self.lengthSquared()


class SpringLink(PhysicalNodeLink):

    def __init__(self,n1:PhysicalNode, n2:PhysicalNode, **kwargs):

        #if self.n1 is None or self.n2 is None:
        super().__init__(n1,n2, **kwargs)

        self.desiredLength = float(kwargs["desiredLength"])
        self.springStrength = float(kwargs["springStrength"])

    def spring_displacement(self, length:float = None):
        if length is None:
            return self.length() - self.desiredLength
        else:
            return self.desiredLength - length

    def get_spring_strength(self, length:float=None):
        return self.spring_displacement(length) * self.springStrength


def polynomialPotential(value, bias:float = 0.0, firstDegreeFactor:float=0.0, secondDegreeFactor:float=0.0, thirdDegreeFactor:float=0.0):
    #todo reavaluate this as is not actually used in my implementation and could be done much more cleanly

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
    nw.gen_nodes_with_radius(50, 20)

    nw.ge_links_between_all(nw, PhysicalNodeLink)

    nw.iterative_calculate_force_vectors(10, 0.0, 0.5)
    nw.iterative_calculate_force_vectors(10, 0.5, 0.0)

    nw.connect_links(nw.gen_closest_n_links(1, nw.nodes, nw.links), SpringLink, desiredLength=1.0, springStrength= 1.0)

    nw.remove_links_of_type(PhysicalNodeLink)

    #nw.connect_node_islands()

    nw.iterative_calculate_force_vectors(10, 0.0, 0.5)
    nw.iterative_calculate_force_vectors(10, 0.5, 0.0)

    nw.getStats()

    Plotter(nw,15)

def plot_quiver():
    print(pplot.get_backend())

    nw = NodeWrangler()

    nw.gen_nodes_with_radius(100, 10)
    nw.ge_links_between_all(1.0)

    nw.iterative_calculate_force_vectors(10, 0.0, 2.0)
    nw.iterative_calculate_force_vectors(10, 2.0, 0.0)


    nw.assign_moveVec_to_displacement_from_origional_position()
    nw.disp_quiver()
    nw.reset_move_vectors()


if __name__ == "__main__":
    test_graph()
