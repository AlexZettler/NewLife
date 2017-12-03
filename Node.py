import math
import random
import copy

import IDtrack
import graphics as gfx


#import matplotlib.pyplot as pplot
#import numpy


class NodeWrangler(object):
    """"
    Generates, tracks and alters all points and links in the graph
    """
    def __init__(self):

        self.nodes = IDtrack.IDTracker()
        self.links = IDtrack.IDTracker()

    def __str__(self):
        return "\n".join(str(self.nodes[n]) for n in self.nodes)

    def get_stats(self) -> str:
        '''
        Returns a bunch of stats regarding the NodeWrangler

        :return: A string containing a bunch of stats of the NodeWrangler
        '''
        return "Getting you some stats of the graph...\n\t" + "\n\t".join(
            ("there are {} links".format(len(self.links)),
             "unlinked nodes: {}".format(self.get_unlinked_nodes()),
             "number of links per node: {}".format(self.get_node_links()),
             "number of islands: {}".format(self.get_node_islands())))

    def generate_nodes(self, n:int, radius:float)->None:
        '''
        Generates n nodes randomly distributed in circle of a certain radius and stores them.

        :param n: Number of nodes to generate
        :param radius: Radius of circle to generate nodes inside of
        :return: None
        '''
        for i in range(n):
            random_node = Pos2D.create_rand_pos_in_circle(radius)
            self.nodes.store_value(PhysicalNode(random_node.x, random_node.y, 1.0))

    '''
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
    '''


    @classmethod
    def generate_links_between_all(self, tracker, link_type, **kwargs)->None:

        #todo rethink this class and implementation. Perhaps divide the generated nodes into sectors and only connect nodes between nodes in current and 8 neighboring sectors, to keep connections local
        # todo This guy can be optimized eventually by making a list that is iterated through and have values poped as the list is iterated through.
        '''
        Generates links between all nodes of the NodeWrangler object. This can be very slow with large numbers of nodes

        :param tracker: The NodeWrangler object who's nodes are to be connected
        :param link_type: The type of link generated
        :param kwargs: Additional arguments to pass into the link generation
        :return: None
        '''

        completedkeys = []

        if isinstance(link_type(Pos2D(0, 0), Pos2D(0, 0), **kwargs), NodeLink):

            for key in tracker.nodes:

                completedkeys.append(key)

                connectToList = {k: tracker.nodes[k] for k in set(tracker.nodes) - set(completedkeys)}

                for connection_key in connectToList:
                    pl = link_type(tracker.nodes[key], tracker.nodes[connection_key], **kwargs)

                    tracker.links.store_value(pl)

                    # print("pl={}".format(pl.get_id()))
                    tracker.nodes[key].links.append(pl.get_id())
                    tracker.nodes[connection_key].links.append(pl.get_id())

        else:
            raise ValueError

    def remove_links_of_type(self, instance_type: type):
        '''
        Iterates through all links and removes any links of a certain type

        :param instance_type: the type of the link that is to be removed
        :return: None
        '''

        #build list of keys to remove
        keysToRemove=[]
        for l in self.links:
            if type(self.links[l]) is instance_type:
                keysToRemove.append(l)
                #print(self.links[selection].n1.links)

        #removes keys
        for l in keysToRemove:
            self.remove_link(l)

    def remove_link(self,link_id)-> None:
        '''
        Removes a link and both nodes that reference it

        :param link_id: id of link to remove
        :return: None
        '''
        self.links[link_id].n1.links.pop(self.links[link_id].n1.links.index(link_id))
        self.links[link_id].n2.links.pop(self.links[link_id].n2.links.index(link_id))
        self.links.delete_value(link_id)

    #Returns a list of all nodes that are not linked
    def get_unlinked_nodes(self):
        return [n for n in self.nodes if len(self.nodes[n].links) == 0]

    # Returns a list of the number of links of each node
    def get_node_links(self):
        return {n: len(self.nodes[n].links) for n in self.nodes}

    def get_node_islands(self) -> {int: []}:
        '''
        Iterates through all nodes and uses a length first search to populate a dictionary containing all node ids

        :return: A dictionary containing a list of nodes for each island key
        '''

        #todo comment this(large task)
        #todo change dictionary to list[islands] of lists[ids]

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

                        c = self.links[l].get_connected(open_set_selection).get_id()

                        if c not in closed_set:

                            open_set.insert(0, c)
                            node_islands[island_counter].append(int(c))
                        closed_set.append(c)

                island_counter += 1


        return node_islands

    def connect_node_islands(self):

        while True:

            ni = self.get_node_islands()

            if len(ni) <= 1:
                break

            island_wrangler = NodeWrangler()

            # create nodegroups with collection of nodes from each island
            for i in ni:
                island_wrangler.nodes.store_value(NodeGroup(self, [n for n in ni[i]]))


            # create a lose link to find shortest connections
            island_wrangler.generate_links_between_all(island_wrangler, NodeLink)


            #connect all closest node island groups
            for island_node in island_wrangler.nodes:
                island_wrangler.connect_links(
                    island_wrangler.gen_closest_n_links(1, island_wrangler.nodes, island_wrangler.links), SpringLink,
                    springStrength=0.0, desiredLength=0.0)

            island_wrangler.remove_links_of_type(NodeLink)

            for island_link in island_wrangler.links:

                first_group = island_wrangler.links[island_link].n1
                second_group = island_wrangler.links[island_link].n2

                NodeGroup.true_closet_connect(self, first_group, second_group)


            island_wrangler.remove_links_of_type(PhysicalNodeLink)

            Plotter(island_wrangler, 5)


    def gen_closest_n_links(self, N, node_tracker, link_tracker):

        retlinks = []

        for n in self.nodes:

            print("n={}/{}total nodes. There are {} links".format(n, len(self.nodes)-1,len(self.nodes[n].links)))

            #creates a dictionary of connections of the current node as LinkKey: LinkLength
            lenDict = {k:self.links[k].length_squared() for k in self.nodes[n].links}

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
        '''
        :param links: a list of links to connect
        :param link_type: the type of link desired
        :param kwargs: additional arguements passed to the constructer of link_type, valid args are as follows:
            PhysicalNodeLink: None
            SpringLink: springStrength:float, desiredLength:float
        :return: None
        '''

        #ensures that the link to create is a nodelink
        if isinstance(link_type(Pos2D(0, 0), Pos2D(0, 0), **kwargs), NodeLink):
            for l in links:

                self.links[l] = link_type(self.links[l].n1, self.links[l].n2, **kwargs)
        else:
            raise TypeError

    '''
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
    '''

    def reset_move_vectors(self):
        for n in self.nodes:
            self.nodes[n].moveVec = Pos2D(0.0, 0.0)

    def calculate_force_vectors_iterative(self,
                                          iterations:int,
                                          start_strength:float,
                                          end_strength:float,
                                          charge_force_multiplier=1.0,
                                          spring_force_multiplier=1.0):
        '''
        Iteratively applies forces to nodes based on connections
        Magnitude of force strength scales linerly from startStr to endStr

        :param iterations:
        :param start_strength: strength of force application at the start of the iteration
        :param end_strength: strength of force application at the end of the iteration
        :return: None
        '''

        #generate
        for x in range(iterations):
            #m is slope
            m = (end_strength - start_strength) / iterations
            #b is offset
            b = start_strength

            mag = m * x + b

            self.calc_force_vectors(mag*charge_force_multiplier, mag*spring_force_multiplier)
            self.apply_vectors()

    def calc_force_vectors(self, charge_strength:float, spring_strength:float):

        '''
        Calculates move vectors of all nodes in a nodewrangler
        Move vectors based on forces applied from links

        :param charge_strength: Value to multiply charge forces by
        :param spring_strength: Value to multiply spring forces by
        :return: None
        '''

        print("Iterating through vectors for force calculation")

        for l in self.links:

            #Hypotenuse stuff
            hypotenuse_squared = self.links[l].length_squared()
            if hypotenuse_squared == 0.0:
                print("n1={}, n2={}".format(self.links[l].n1, self.links[l].n2))
            hyp = math.sqrt(hypotenuse_squared)


            #ang is the degree measure between the displacement vecotr of the link and the positive x axis
            ang = math.atan2(self.links[l].delta_y, self.links[l].delta_x)


            #How much force is applied to the object
            move_magnatude = 0.0

            # need to recalculate hypsquared here so optional parameter to ease calculation
            if isinstance(self.links[l], PhysicalNodeLink):
                move_magnatude -= charge_strength * self.links[l].charged_particle_force(hypotenuse_squared)

            if isinstance(self.links[l],SpringLink):
                #print("spring strength applied")
                move_magnatude += spring_strength * self.links[l].get_spring_strength(hyp)


            print("the movement vector for link: {}"
                  "is in direction: {}"
                  "and has a magnatude of: {}".format(l,to_degree(ang),move_magnatude))

            #print("cosang: {} sinang: {} chargeMul: {} ang: {}".format(coscomp*180/math.pi,sincomp*180/math.pi,chargemul,ang*180/math.pi))
            #print("N1:{},N2:{} ".format(self.links[l].n1.moveVec,self.links[l].n2.moveVec))

            self.links[l].n1.move_vector.x += -move_magnatude / self.links[l].n1.weight * math.cos(ang)
            self.links[l].n1.move_vector.y += -move_magnatude / self.links[l].n1.weight * math.sin(ang)

            self.links[l].n2.move_vector.x += move_magnatude / self.links[l].n2.weight * math.cos(ang)
            self.links[l].n2.move_vector.y += move_magnatude / self.links[l].n2.weight * math.sin(ang)

    def add_random_force_vectors(self,radius):

        for n in self.nodes:
            if isinstance(n,PhysicalNode):
               n.move_vector += Pos2D.create_rand_pos_in_circle(radius)

    def apply_vectors(self):

        #Apply the moveVectors
        for n in self.nodes:

            if self.nodes[n].rooted == False:

                #print("type: {}, node x: {}, y: {}".format(type(self.nodes[n]),self.nodes[n].x, self.nodes[n].y)) #+= self.nodes[n].moveVec.x
                #print("type: {}, movevec x: {}, y: {}".format(type(self.nodes[n].moveVec),self.nodes[n].moveVec.x, self.nodes[n].moveVec.y))#+= self.nodes[n].moveVec.y
                self.nodes[n].x += self.nodes[n].move_vector.x
                self.nodes[n].y += self.nodes[n].move_vector.y

        self.reset_move_vectors()

    '''
    def disp_plot(self):
        data = self.gen_quiver_graph_data()
        pplot.scatter(data[0], data[1])
        pplot.show()

    def disp_quiver(self):
        data = self.gen_quiver_graph_data()
        pplot.quiver(data[0], data[1], data[2], data[3], label = "This label")
        pplot.show()

    def gen_quiver_graph_data(self):
        xdata = [self.nodes[n].x for n in self.nodes]
        ydata = [self.nodes[n].y for n in self.nodes]

        xmoveData = [self.nodes[n].moveVec.x for n in self.nodes]
        ymoveData = [self.nodes[n].moveVec.y for n in self.nodes]

        data = numpy.array([xdata, ydata,xmoveData,ymoveData])

        print("some graph data:\n\t{},{} pieces of data in xdata,ydata\n\t{},{} pieces of data in xmovedata,ymovedata".format(len(xdata),len(ydata),len(xmoveData),len(ymoveData)))

        #p = pplot.quiver(data[0], data[1], data[2], data[3])
        return data
    '''

class Plotter(object):

    '''
    An object used to draw a Nodewrangler object into a graphics window
    '''

    def __init__(self, node_wrangler_object: NodeWrangler, scale):
        '''
        :param node_wrangler_object: The nodewrangler object to be plotted
        :param scale: How distant the points should be drawn
        '''

        self.canvas = gfx.GraphWin("WorldGen", 1000, 800)

        if isinstance(node_wrangler_object.nodes[0], PhysicalNode):
            nodestoplot = node_wrangler_object.nodes

        else:
            raise TypeError

        self.plot_points(nodestoplot, scale)

        #print(nwo.nodes.keys())

        linestoplot = {k: [v.n1, v.n2] for k, v in node_wrangler_object.links.items()}
        self.plot_lines(linestoplot, scale)

        self.canvas.getMouse()

    def plot_points(self, point_id_dict: dict, scale=1.0):
        # todo commecnt this
        '''

        :param point_id_dict:
        :param scale:
        :return:
        '''

        for n in point_id_dict:

            cur_node = point_id_dict[n]
            c = gfx.Circle(center=self._generate_point(cur_node.x, cur_node.y, scale),radius=scale*cur_node.weight*0.5)

            c.draw(self.canvas)


    def _generate_point(self, x, y, scale) -> gfx.Point:
        '''
        :param x: x position of the point
        :param y: y position of the point
        :param scale: How much to zoom 
        :return:
        '''

        return gfx.Point(scale * x + self.canvas.width // 2, scale * y + self.canvas.width // 2)

    def plot_lines(self, line_id_dict: dict, scale=1.0) -> None:
        '''
        :param line_id_dict:
        :param scale: How 
        :return: None
        '''

        for l in line_id_dict:
            #print(point_id_dict[0][0])
            cur_line = line_id_dict[l]
            p1 = cur_line[0]
            p2 = cur_line[1]

            line = gfx.Line(self._generate_point(p1.x, p1.y, scale), self._generate_point(p2.x, p2.y, scale))

            line.draw(self.canvas)


#Pos2D class Collection
class Pos2D(object):

    '''
    Pos2D:
        represents a x,y pair of co-ordinates with override methods
    '''

    def __init__(self, x, y):
        print("position 2D")
        self.x, self.y = float(x), float(y)

    def __str__(self):
        return "Pos2D-x:{}, y:{}".format(self.x, self.y)


    '''

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

    def __sub__(self, other):
        #todo perhaps implement some sort of child member add functionality? Is this a good idea?
        return Pos2D(self.x - other.x, self.y - other.y)

    '''

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
    def create_rand_pos_in_circle(cls, radius=1.0):

        '''
        return a point randomly generated that is within the circle of given radius

        :param radius: the radius of the circle to generate the point within
        :return: A Pos2D object with it's coordinates inside the circle of given radius
        '''

        point_outside_circle = True
        # Continue loop until loop is inside circle
        while point_outside_circle:
            x = (random.random() - 0.5) * radius
            y = (random.random() - 0.5) * radius

            #checks if point is inside circle
            if x ** 2 + y ** 2 < radius ** 2:
                point_outside_circle = False

        return Pos2D(x, y)

class Node(Pos2D, IDtrack.TrackedObject):
    '''
    A tracked position with a list of connection ids
    '''
    def __init__(self, x:float, y:float):
        '''
        :param x: x position of the Node
        :param y: y position of the Node
        '''

        #Explicitly call Pos2D's init method because of multiple inheritances
        Pos2D.__init__(self,x,y)

        #A list of link ids that contain this Node
        self.links = []

    def track(self, id):
        IDtrack.TrackedObject.__init__(self, id)


    def __str__(self):
        if id is None:
            return super().__str__()
        else:
            return "id: {}, x: {}, y: {}".format(self.get_id(), self.x, self.y)

    def __int__(self):
        return self.get_id()

class PhysicalNode(Node):

    '''
    A Node that has physics attributes

    :attribute move_vector: Pos2D component representing forces applied based on connections every tick
    :attribute weight: A float representing charge and inertia
    :attribute rooted: A boolean representing if the PhysicalNode object will move when forces are applied to it
    '''

    def __init__(self, x:float, y:float, weight:float=1.0):
        Node.__init__(self, x, y)

        self.move_vector = Pos2D(0.0, 0.0)
        self.weight = weight
        self.rooted = False

    def __str__(self):
        return "{}, moveVec: {}".format(str(super), self.move_vector)

    def root(self,state:bool):
        self.rooted = state


class NodeGroup(PhysicalNode):

    '''
    NodeGroup represents a collection of nodes with a position at the average position of the nodes

    :attribute contained_nodes: A list of nodes in the parent NodeWrangler
    '''

    def __init__(self, node_wrangler_parent: NodeWrangler, nodes: list):
        '''
        :param node_wrangler_parent: the main NodeWrangler object that the list of nodes are tracked from
        :param nodes: The nodes that are part of the group
        '''

        self.contained_nodes = nodes
        center = Pos2D(0.0, 0.0)

        for k in self.contained_nodes:
            center += node_wrangler_parent.nodes[k]
        # print(len(Nodes))
        center /= float(len(nodes))

        super().__init__(center.x, center.y)

    @classmethod
    def true_closet_connect(self, node_wrangler_parent: NodeWrangler, first_group, second_group)->None:
        '''
        :param node_wrangler_parent: 
        :param first_group: A NodeGroup object representing a collection of nodes
        :param second_group: A NodeGroup object representing a collection of nodes
        :return: None
        '''

        shortest_link = None

        check_link = None

        print(first_group.contained_nodes)
        print(second_group.contained_nodes)


        #iterate through all possible node connections
        for n1 in first_group.contained_nodes:

            for n2 in second_group.contained_nodes:

                # print("n1={},n2={}".format(n1, n2))

                check_link = NodeLink(node_wrangler_parent.nodes[n2], node_wrangler_parent.nodes[n1])

                if shortest_link is None:
                    shortest_link = check_link

                elif check_link.length_squared() < shortest_link.length_squared():

                    shortest_link = check_link

        if check_link is not None:

            check_link = SpringLink(check_link.n1, check_link.n2, desiredLength=1.0, springStrength= 1.0)

            node_wrangler_parent.links.store_value(check_link)
            node_wrangler_parent.nodes[check_link.n1.get_id()].links.append(check_link.get_id())
            node_wrangler_parent.nodes[check_link.n2.get_id()].links.append(check_link.get_id())

            first_group.contained_nodes += second_group.contained_nodes

        else:
            raise ValueError


#NodeLink class Collection
class NodeLink(IDtrack.TrackedObject):
    '''
    Represents a pair of (bidirectional) connected nodes

    n1 is first node
    n2 is second node
    '''

    def __init__(self,n1,n2,**kwargs):
        self.n1, self.n2 = n1, n2
        if n1 == n2:
            raise ValueError

        #the init call must be called here
        IDtrack.TrackedObject.__init__(self, None)

    def length_squared(self):
        return self.delta_x ** 2 + self.delta_y ** 2

    def length(self):
        return math.sqrt(self.length_squared())

    def __len__(self):
        return self.length()

    @property
    def delta_x(self):
        return self.n2.x - self.n1.x

    @property
    def delta_y(self):
        return self.n2.y - self.n1.y


    def get_line_equation(self):
        '''
        calculates the linear equation of the link line

        :return: a tupple containing the m,b pair.
        '''

        m = self.delta_x/self.delta_y
        b = self.n1.y - m*self.n1.x

        #todo return the range of x,y values for the given line

        return (m, b)

    def get_connected(self, n):
        '''
        Returns the id of the connected node given the id of a node in a link

        :param n: id of link you wish to find connection of
        :return: id of connected node
        '''
        intn1 = int(self.n1)
        intn2 = int(self.n2)
        intn = int(n)

        if intn == intn1:
            retid = self.n2
        elif intn == intn2:
            retid = self.n1
        else:
            raise KeyError("The two nodes requesting connection have the same ID")

        #print("Printing the connected link of node {}\n\tThe connected node is: {}\n\tThe link id is: {}".format(self.n1, self.n2, retid))
        return retid


class PhysicalNodeLink(NodeLink):
    def __init__(self,n1, n2, **kwargs):

        super().__init__(n1, n2, **kwargs)

    def charged_particle_force(self, lenSquared:float = None):

        if lenSquared is not None:
            return self.n1.weight * self.n2.weight / lenSquared
        else:
            return self.n1.weight * self.n2.weight / self.length_squared()


class SpringLink(PhysicalNodeLink):

    def __init__(self,n1:PhysicalNode, n2:PhysicalNode, **kwargs):

        #if self.n1 is None or self.n2 is None:
        super().__init__(n1,n2, **kwargs)


        #desiredLength is the distance that the springlink rests at
        self.desiredLength = float(kwargs["desiredLength"])

        #springStrength is how much force is applied per unit from the spring's rest length
        self.springStrength = float(kwargs["springStrength"])

    def spring_displacement(self, length:float = None):
        if length is None:
            return self.desiredLength - self.length()
        else:
            return self.desiredLength - length

    def get_spring_strength(self, length:float=None):
        return self.spring_displacement(length) * self.springStrength

def to_degree(radian):
    return radian*180/math.pi

def to_radian(degree):
    return degree / 180*math.pi

def test_graph():

    nw = NodeWrangler()
    nw.generate_nodes(20, 100)

    Plotter(nw, 5.0)

    nw.generate_links_between_all(nw, PhysicalNodeLink)

    #nw.iterative_calculate_force_vectors(10, 0.0, 0.5)
    #nw.iterative_calculate_force_vectors(10, 0.5, 0.0)

    nw.connect_links(nw.gen_closest_n_links(1, nw.nodes, nw.links), SpringLink, desiredLength=1.0, springStrength= 1.0)

    nw.remove_links_of_type(PhysicalNodeLink)

    nw.connect_node_islands()

    nw.calculate_force_vectors_iterative(10, 0.0, 0.5, charge_force_multiplier=1.0, spring_force_multiplier=1.0)
    nw.calculate_force_vectors_iterative(10, 0.5, 0.0, charge_force_multiplier=1.0, spring_force_multiplier=1.0)

    print(nw.get_stats())
    #Plotter(nw,15)

    #nw.calc_force_vectors(1.0,0.2)
    nw.apply_vectors()

    #Plotter(nw, 15)

    #nw.calculate_force_vectors_iterative(10, 0.5, 0.2)

    #nw.add_random_force_vectors(5.0)
    nw.apply_vectors()

    Plotter(nw,5.0)

'''
def plot_quiver():
    print(pplot.get_backend())

    nw = NodeWrangler()

    nw.generate_nodes(100, 10)
    nw.gen_links_between_all(nw, NodeLink)

    nw.iterative_calculate_force_vectors(10, 0.0, 2.0)
    nw.iterative_calculate_force_vectors(10, 2.0, 0.0)


    nw.assign_moveVec_to_displacement_from_origional_position()
    nw.disp_quiver()
    nw.reset_move_vectors()
'''

if __name__ == "__main__":
    test_graph()
