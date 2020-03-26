from src.strand import bond
from src.util import util
import queue
import copy
from collections import defaultdict


class BondGraph:
    """
    Data structure of the bond graph, undirected graph with loops
    """

    V = []
    '''vertices'''

    E = []
    '''edges'''

    adj = []
    '''adjacency list'''

    hidden = []
    '''hidden domains of strands'''

    loop = []
    '''loops in the graph'''

    species = []
    '''species of the graph'''

    speciesnum = 0
    '''species number'''

    color = []
    '''inherit from StrandGraph'''

    def __init__(self, V, color, E):
        self.V = V
        self.color = color
        self.E = E
        self.adj = [[] for _ in self.V]
        self.hidden = []
        self.loop = []
        self.species = []

    def create_bond(self, v1, v2, d1, d2):
        """
        create bond

        :param v1: the strand to be examined
        :param v2: the target strand to check
        :param d1: the domain position on the examining strand of the bond
        :param d2: the domain position on the target strand of the bond
        :return:
        """
        flag = 0
        if len(self.adj[v1]) != 0:
            for i in self.adj[v1]:
                if v2 == i.node2:
                    flag = 1
                    break
            # if there are bonds between two strands
            if flag == 1:
                i.appenddom(d1, d2)
            # new bond
            else:
                self.adj[v1].append(bond.Bond(v1, v2, [d1], [d2]))
        # new bond
        else:
            self.adj[v1].append(bond.Bond(v1, v2, [d1], [d2]))

    def add_edges(self, v1, v2, d1, d2):
        """
        add edges to the bond graph

        :param v1: strand 1
        :param v2: strand 2
        :param d1: domain 1, corresponds to strand 1
        :param d2: domain 2, corresponds to strand 2
        """
        self.create_bond(v1, v2, d1, d2)
        self.create_bond(v2, v1, d2, d1)
        #self.map_colors({(v1, d1), (v2, d2)})

    def check_in_loop(self, startnode, endnode):
        """
        check if the two vertices are in the same loop

        :param startnode: start vertex
        :param endnode: end vertex
        :return: True if in the same loop, False otherwise
        """
        for i in self.loop:
            if startnode in i and endnode in i:
                return True
        return False

    def check_junction(self, startnode, endnode):
        """
        check if the bond between startnode and endnode is a part of junction

        :param startnode:
        :param endnode:
        :return:
        """
        for i in self.loop:
            if startnode in i and endnode in i:
                for j in range(len(i)):
                    if j == 0:
                        pre = len(i) - 1
                        suc = j + 1
                    elif j == len(i) - 1:
                        pre = j - 1
                        suc = 0
                    else:
                        pre = j - 1
                        suc = j + 1
                    n = []

                    for b in self.adj[i[j]]:
                        if b.node2 == i[pre] or b.node2 == i[suc]:
                            n.append(b.dom)
                    if not util.check_continuity(n[0], n[1]):
                        return False
                return True
        return False

    def get_bond(self, startstrand, endstrand):
        """
        get the bond between the two vertices

        :param startstrand: start strand
        :param endstrand: end strand
        :return: get the bond between two strands
        """
        for i in self.adj[startstrand]:
            if i.node2 == endstrand:
                return i

    def search_path(self, startstrand, stack, endstrand, edge):
        """
        search a path from the start node to the end node in the spanning tree
        """
        stack.append(startstrand)

        if endstrand == stack[len(stack) - 1]:
            self.loop.append(stack)
            return True
        else:
            for bond in self.adj[startstrand]:
                if bond in edge and bond.node2 not in stack:
                    if self.search_path(bond.node2, stack, endstrand, edge):
                        return True
            stack.pop()
            return False

    def spanning(self, mark, edge, startstrand, times, depth):
        """
        building a spanning tree of the bond graph

        :param mark: list of boolean numbers of vertices to show if they have been visited
        :param edge: edges of a spanning tree of the bond graph
        :param startstrand: starting strand to find spanning tree
        :param times: list of integers of vertices indicating the times of this function is called
        :param depth: current cycle number of calling this function
        """
        for bond in self.adj[startstrand]:
            if mark[bond.node2] is False:
                edge.append(bond)
                for b in self.adj[bond.node2]:
                    if b.node2 == bond.node1:
                        edge.append(b)
                        continue
                mark[bond.node2] = True
                times[bond.node2] = depth
                self.spanning(mark, edge, bond.node2, times, depth)

    def call_spanning(self, mark, edge, times, depth):
        """
        find a spanning forest of the bond graph and store the loops existing in it

        :param mark: list of boolean numbers of vertices to show if they have been visited
        :param edge: edges of a spanning tree of the bond graph
        :param times: list of integers of vertices indicating the times of spanning() is called
        :param depth: current cycle number of calling spanning()
        """
        for i in range(0, len(self.V)):
            if not mark[i]:
                mark[i] = True
                times[i] = depth
                self.spanning(mark, edge, i, times, depth)

                for j in range(0, len(self.V)):
                    if mark[j] and times[j] == depth:
                        for bond in self.adj[j]:
                            if mark[bond.node2]:
                                if bond not in edge:
                                    # there is a loop
                                    # no need to store the loop in cycle if the loop is already found
                                    if self.check_in_loop(j, bond.node2):
                                        continue
                                    self.search_path(j, [], bond.node2, edge)
                depth += 1
        self.species = times
        self.speciesnum = depth

    def find_loops(self):
        """
        find loops in the bond graph
        """
        mark = [False for _ in self.V]
        edge = []
        times = [0 for _ in self.V]
        # find loops in bond graph and store them in cycle[]
        self.call_spanning(mark, edge, times, 0)

    def store_hidden(self):
        """
        store all the hidden domains in hidden[]
        """
        for i in self.adj:
            # the case of harpin loop and the case multiple bonds exist between two strands
            for b in i:
                if b.node2 == b.node1 or len(b.dom) > 1:
                    minn = min(b.dom)
                    maxn = max(b.dom)
                    for j in range(minn + 1, maxn):
                        self.hidden.append((b.node1, j))

        for j in self.loop:
            loopedge = []
            for k in range(0, len(j)):
                if k == len(j) - 1:
                    loopedge.append(self.get_bond(j[k], j[0]))
                    loopedge.append(self.get_bond(j[k], j[k - 1]))
                elif k == 0:
                    loopedge.append(self.get_bond(j[k], j[k + 1]))
                    loopedge.append(self.get_bond(j[k], j[len(j) - 1]))
                else:
                    loopedge.append(self.get_bond(j[k], j[k + 1]))
                    loopedge.append(self.get_bond(j[k], j[k - 1]))
            for p in range(0, len(j)):
                minn = max(loopedge[p * 2].dom)
                maxn = max(loopedge[p * 2 + 1].dom)
                if minn > maxn:
                    b = minn
                    minn = maxn
                    maxn = b

                for x in range(minn + 1, maxn):
                    if (loopedge[p * 2].node1, x) not in self.hidden:
                        self.hidden.append((loopedge[p * 2].node1, x))

    def check_bonded(self, node):
        """
        check if the node is bonded to another node

        :param node: a domain on a strand
        :return: the domain it binds to if True, otherwise None
        """
        v, n = node

        for i in self.adj[v]:
            for j in range(0, len(i.dom)):
                if n == i.dom[j]:
                    return i.node2, i.dom2[j]
        return None

    def check_strand_is_bonded(self, v):
        """
        check if the strand is bonded

        :param v: strand
        :return: domains on strand that is bonded
        """
        bonddomains = []

        for i in self.adj[v]:
            bonddomains += i.dom
        return bonddomains

    def get_species(self):
        """
        get all species in one bondgraph

        :return:all species
        """
        speciesnodes = [[] for _ in range(0, self.speciesnum)]
        for j in range(0, len(self.V)):
            speciesnodes[self.species[j]].append(self.V[j])
        return speciesnodes

    def get_connection(self, startnode, endstrand):
        """
        get connection nodes from the startnode to the endstrand

        :param startnode: starting node
        :param endstrand: end strand
        :return: path from startnode to end strand
        """
        q = queue.Queue()
        visited = [False for _ in self.V]
        flag = False
        prev = {}
        connect = []

        q.put(startnode)

        while not q.empty():
            curstrand, curdomain = q.get()
            visited[curstrand] = True

            if curstrand == endstrand:
                flag = True
                break

            for i in self.adj[curstrand]:
                if not visited[i.node2]:
                    d2 = max(i.dom2)
                    q.put((i.node2, d2))
                    prev[(i.node2, d2)] = (curstrand, curdomain)

        if not flag:
            return None

        key = (curstrand, curdomain)
        while key in prev:
            connect.append(key)
            connect.append(prev[key])
            key = prev[key]

        return connect

    def get_direction(self, startstrand, endstrand):
        """
        Note that we don't consider loops here.
        :param endstrand:
        :param startstrand:
        :return:
        """

        if len(self.V) <= 2:
            return False

        q = queue.Queue()
        visited = [False for _ in range(len(self.V))]
        direction = [0 for _ in range(len(self.V))]
        cnt = 0

        q.put(startstrand)
        while not q.empty():
            cnt += 1
            curstrand = q.get()

            if curstrand == endstrand:
                break

            # Exception occurs
            if len(self.adj[curstrand]) == 0:
                continue

            for i in self.adj[curstrand]:
                if not visited[i.node2]:
                    q.put(i.node2)
                    direction[i.node2] = util.flip(direction[i.node1])

        if direction[startstrand] != direction[endstrand] and cnt > 2:
            return True
        return False

    def merge_bonds_ignoring_nodes(self, v, nodes):
        """
        merge bonds on v ignoring bonds connected to nodes

        :param v: a strand
        :param nodes: a set of strands
        :return: bonds
        """
        if not isinstance(v, tuple):
            return self.adj[v]
        bonds = []
        for i in v:
            for j in self.adj[i]:
                if j.node2 not in nodes:
                    bonds.append(j)
        return bonds


class SubBondGraph(BondGraph):
    """
    Subgraph of bondgraph
    """

    Vsub = []
    '''vertices'''

    adj = []
    '''adjacency list'''

    colormap = defaultdict(set)
    '''map colors to vertices'''

    colorset = set()
    '''set of colors'''

    def __init__(self, V, color, adj, Vp):
        self.Vsub = V
        self.color = color
        self.colorset = set()
        self.colormap = defaultdict(set)

        BondGraph.__init__(self, Vp, color, [])
        if len(self.Vsub) == 1:
            self.colormap[self.color[self.Vsub[0]]].add(self.Vsub[0])
            self.colorset.add(self.color[self.Vsub[0]])
            return
        self.adj = copy.copy(adj)

        visited = [False for _ in V]
        d = dict(zip(V, visited))

        for v in Vp:
            if v not in V:
                self.adj[v] = []
                continue
            for i in range(0, len(self.adj[v])):
                if self.adj[v][i].node2 not in V:
                    self.adj[v].pop(i)
                else:
                    if d[self.adj[v][i].node2]:
                        continue
                    else:
                        self.map_colors({(v, self.adj[v][i].dom[0]), (self.adj[v][i].node2, self.adj[v][i].dom2[0])})
                        d[v] = True

    def map_colors(self, edge):
        for x in edge:
            self.colormap[self.color[x[0]]].add(x[0])
            self.colorset.add(self.color[x[0]])
