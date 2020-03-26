from src.util import util
from src.util import cexception as ex
from src.strand import bond_graph
import copy


class StrandGraph:
    """
    Data structure of the strand graph
    """

    V = []
    '''vertices'''

    length = []
    '''length of the vertices (in length of domains)'''

    color = []
    '''strand type'''

    A = []
    '''admissible edges, edges that bonds can exist in'''

    toehold = []
    '''a list of dictionaries to store if toehold exist in an edge'''

    E = []
    '''actual edges, edges that already exist'''

    strands = []
    '''strands in this graph'''

    bondgraph = None
    '''bond graph derived from this graph'''

    def __init__(self, strands):
        """
        build strand graph for strands

        :param strands: list of Strand type objects
        """
        def comp(a, b):
            if a[2] == b[2]:
                return False
            else:
                return True

        V = [i for i in range(0, len(strands))]

        color = [strands[i].color for i in range(0, len(strands))]
        length = [len(strands[i].domains) for i in range(0, len(strands))]
        A = []
        E = []
        toehold = {}

        tableS = []
        i = 0
        for s in strands:
            j = 0
            for d in s.domains:
                tableS.append([(d.name, (i, j), d.comp, d.bond, d.bondname, d.toehold), False])
                j += 1
            i += 1

        for i in range(0, len(tableS)):
            cur = tableS[i]
            if cur[1]:
                bondexist = True
            else:
                bondexist = False

            for j in range(i + 1, len(tableS)):
                if tableS[j][0][0] == cur[0][0] and comp(tableS[j][0], cur[0]):
                    A.append({cur[0][1], tableS[j][0][1]})
                    if tableS[j][0][3] == cur[0][3] == 1 and tableS[j][0][4] == cur[0][4] and not tableS[j][1] and not cur[1]:
                        E.append({cur[0][1], tableS[j][0][1]})
                        bondexist = True
                        tableS[j][1] = True
                        cur[1] = True
                    if tableS[j][0][5] == cur[0][5] == True:
                        toehold[frozenset({cur[0][1], tableS[j][0][1]})] = True
                        # toehold.append({'e': {cur[0][1], tableS[j][0][1]}, 't': True})
                    else:
                        toehold[frozenset({cur[0][1], tableS[j][0][1]})] = False
                        # toehold.append({'e': {cur[0][1], tableS[j][0][1]}, 't': False})
            if not bondexist and cur[0][3]:
                raise ex.SpeciesError("species text representation error in domain " + cur[0][0])

        self.V = V
        self.length = length
        self.color = color
        self.A = A
        self.toehold = toehold
        self.E = E
        self.strands = strands
        self.build_bond_graph()

    def reconstruct(self, E):
        """
        reconstruct the graph using the set of edges

        :param E: edges
        """
        self.E = E
        self.build_bond_graph()

    def delete_vertex(self, v):
        strands = copy.copy(self.strands)
        strands.pop(v)

        self.__init__(strands)

    def get_domain_obj(self, domain):
        """

        :param domain:
        :return:
        """
        return self.strands[domain[0]].domains[domain[1]]

    def check_complementary(self, domain1, domain2):
        """

        :param domain1:
        :param domain2:
        :return:
        """
        obj1 = self.get_domain_obj(domain1)
        obj2 = self.get_domain_obj(domain2)
        if obj1.name == obj2.name:
            if obj1.comp != obj2.comp:
                return True
        return False

    def anchored(self, e):
        """
        check if e is anchored

        :param e: an edge
        :return: True if anchored, False otherwise
        """
        # TODO: CHECK DEFINITION OF ANCHOR
        v, n = util.get_edge_info(e)

        for i in self.bondgraph.adj[v[0]]:
            if n[0] + 1 in i.dom or n[0] - 1 in i.dom:
                if i.node2 == v[1]:
                    return True

        if self.bondgraph.check_junction(v[0], v[1]):
            return True

        return False

    def build_bond_graph(self):
        """
        build bond graph for the strand graph
        """
        graph = bond_graph.BondGraph(self.V, self.color, self.E)
        v = [0 for i in range(2)]
        n = [0 for i in range(2)]

        # insert test data
        # E = [{(0, 3), (1, 1)}, {(0, 4), (2, 1)}, {(0, 2), (3, 3)}, {(1, 2), (3, 0)}, {(0, 1), (3, 5)}, {(4, 1), (4, 3)}]
        for i in self.E:
            t = 0
            for x in i:
                v[t], n[t] = x
                t += 1
            graph.add_edges(v[0], v[1], n[0], n[1])

        graph.find_loops()
        graph.store_hidden()

        self.bondgraph = graph

    def anti_parallel(self, node1, node2):
        """
        check if the bond satisfies the anti-parallel schema

        :param node1:
        :param node2:
        :return: True if satisfy, False otherwise
        """
        node1bdomains = self.bondgraph.check_strand_is_bonded(node1[0])
        if node1bdomains is None:
            return True

        for i in node1bdomains:
            connect = self.bondgraph.get_connection((node1[0], i), node2[0])
            if connect is None:
                continue
            elif len(connect) == 0:
                continue
            else:
                if ((node1[1] > i) and (node2[1] < connect[0][1])) \
                        or ((node1[1] < i) and (node2[1] > connect[0][1])):
                    return True
                else:
                    node2bdomains = self.bondgraph.check_strand_is_bonded(node2[0])
                    if util.get_free_domains([node1[1], i], node1bdomains, node1[1]) > 1 \
                            or util.get_free_domains([node2[1], connect[0][1]], node2bdomains, node2[1]) > 1:
                        return True
                    elif self.bondgraph.get_direction(node1[0], node2[0]):
                        return True
                    return False

        return True

    def available(self, e):
        """
        check if a binding is available for binding

        :param e: a bond
        :return: True if available, False otherwise
        """
        if self.hidden(e):
            return False

        for x in e:
            if self.bondgraph.check_bonded(x) is not None:
                return False

        v, n = util.get_edge_info(e)

        if self.anti_parallel((v[0], n[0]), (v[1], n[1])):
            return True
        '''
        # check if the corresponding two strands are already bonded
        existing0 = []
        existing1 = []
        for bond in self.bondgraph.adj[v[0]]:
            if bond.node2 == v[1]:
                existing0 = bond.dom
                existing1 = bond.dom2

        if len(existing0) == 0:
            if self.check_adjacent_bond(v, n):
                return True
            return False

        # if an admissible egde can exist between the bonds of two strands
        if (n[0] > max(existing0) and n[1] < min(existing1)) or (n[0] < min(existing0) and n[1] > max(existing1)) \
                or (max(existing0) > n[0] > min(existing0) and max(existing1) > n[1] > min(existing1)):
            return True
        else:
            if len(existing0) > 1:
                for i in range(0, len(existing0) - 1):
                    if existing0[i] + 1 in existing0 or existing0[i] - 1 in existing0:
                        continue
                    else:
                        return True
                for i in range(0, len(existing1) - 1):
                    if existing1[i] + 1 in existing1 or existing1[i] - 1 in existing1:
                        continue
                    else:
                        return True
        '''
        return False

    def hidden(self, e):
        """
        check if domains in the bond are hidden in a cycle of strands

        :param e: bond
        :return: True if hidden, False otherwise
        """
        node = []
        flag = 0
        for x in e:
            if x in self.bondgraph.hidden:
                flag += 1
                node.append(x[0])
        if 1 >= flag > 0:
            return True
        elif flag == 2:
            # if both of the domains are hidden and in the same loop, then they are available to bind
            if self.bondgraph.check_in_loop(node[0], node[1]):
                return False
            return True
        else:
            return False

    def check_toehold(self, e):
        """
        check if the binding domains are toehold domains

        :param e: a bond
        :return: True if are toehold domains, False otherwise
        """
        if self.toehold[frozenset(e)]:
            return True
        return False

    def same_species(self, e):
        """
        check if the edge is from the same species

        :param e: edge
        :return: True if it is, False otherwise
        """
        v, n = util.get_edge_info(e)
        if self.bondgraph.species[v[0]] == self.bondgraph.species[v[1]]:
            return True
        return False

    def check_bonded(self, strand1, strand2):
        """
        check if strand1 is bonded to strand2

        :param strand1:
        :param strand2:
        :return: bonded domains in format [{(strand1, dom1), (strand2, dom2)}, ...]
        """
        bondeddomain = []
        for i in self.bondgraph.adj[strand1]:
            if i.node2 == strand2:
                for j in range(0, len(i.dom)):
                    bondeddomain.append({(strand1, i.dom[j]), (strand2, i.dom2[j])})
        return bondeddomain

    def get_domain_comp(self, domain):
        return self.strands[domain[0]].domains[domain[1]].comp

    def check_candidate(self, e1):
        """
        check if e1 is a candidate for rule migration

        :param e1: an admissible edge
        :return:
        """
        notbond = None
        r = []
        rp = []
        tbr = []

        for x in e1:
            bondto = self.bondgraph.check_bonded(x)
            if bondto is None:
                notbond = x
            else:
                tbr.append({x, bondto})
                rp.append(bondto)

        if len(tbr) == 2:
            r = rp
        if len(tbr) == 1:
            r = rp[0]

        return notbond, tbr, r

    def check_switchable(self, notbond, tbr, i):
        """
        check if it is possible to do 3-way migration

        :param notbond: notbond domain
        :param r: replace domain
        :return: True if possible, False otherwise
        """
        potbondto = list(tbr & i)
        potbondto = potbondto[0]

        if notbond in self.bondgraph.hidden:
            return False

        prevE = copy.copy(self.E)
        self.delete_edge(tbr, prevE)
        connect = self.bondgraph.get_connection(notbond, potbondto[0])

        if connect is None:
            self.reconstruct(prevE)
            return False

        self.reconstruct(prevE)

        nbbondto = []
        for i in self.bondgraph.adj[notbond[0]]:
            for j in range(0, len(i.dom)):
                nbbondto.append({(notbond[0], i.dom[j]), (i.node2, i.dom2[j])})

        for i in nbbondto:
            flag = False
            if self.check_toehold(i):
                for j in i:
                    if j[0] != notbond[0]:
                        ip = copy.copy(i)
                        ip.remove(j)
                        td = list(ip)
                        td = td[0]

                        if len(connect) != 0:
                            if (notbond[1] > td[1] and potbondto[1] < connect[0][1]) \
                                    or (notbond[1] < td[1] and potbondto[1] > connect[0][1]):
                                flag = True

                        else:
                            if (notbond[1] > td[1] and potbondto[1] < j[1]) \
                                    or (notbond[1] < td[1] and potbondto[1] > j[1]):
                                flag = True

                        if abs(notbond[1] - j[1]) > 1:
                            flag = True

                    if flag:
                        return True
                    else:
                        continue

        return False

    def delete_edge(self, e, prevE):
        # TODO: CHECK
        E = copy.copy(prevE)
        if e in E:
            E.remove(e)
        self.reconstruct(E)

    def delete_edges_regarding_v(self, v, prevE):
        """
        delete a vertex v and all edges link to it in strandgraph

        :param v:
        :param prevE:
        """
        E = copy.copy(prevE)
        t = 0
        while 1:
            for j in E[t]:
                if v == j[0]:
                    E.pop(t)
                    t -= 1
                    break
            t += 1
            if t >= len(E):
                break
        self.reconstruct(E)

    def check_switchable_2(self, r1, tbr):
        """
        check if it is possible for 4-way migration

        :param r1: replacing bond
        :param tbr: bond to be replaced
        :return: True if possible, False otherwise
        """
        cdomain = list(tbr[0] - r1)
        cdomain = cdomain[0]

        bdomain = list(r1 & tbr[0])
        bdomain = bdomain[0]

        prevE = copy.copy(self.E)

        self.delete_edges_regarding_v(cdomain[0], prevE)

        if self.bondgraph.speciesnum > 2:
            self.reconstruct(prevE)
            return False

        else:
            if self.check_switchable(bdomain, tbr[1], r1):
                self.reconstruct(prevE)
                return True

        self.reconstruct(prevE)
        return False

    def check_hidden_all(self, edge):
        """

        :param edge:
        :return:
        """
        for e in edge:
            v, n = util.get_edge_info(e)
            if (v[0], n[0]) in self.bondgraph.hidden \
                    or (v[1], n[1]) in self.bondgraph.hidden:
                continue
            else:
                return False
        return True

    def get_connect_toehold(self, edge, notbond, potbondto):
        """
        get the toeholds that connect the two domains of edge

        :param edge:
        :param notbond:
        :param potbondto:
        :return:
        """
        v, n = util.get_edge_info(edge)
        strand1 = v[0]
        strand2 = v[1]

        pottoehold = []
        if strand1 != strand2:
            for b in self.bondgraph.adj[strand1]:
                for i in range(0, len(b.dom)):
                    e = {(strand1, b.dom[i]), (b.node2, b.dom2[i])}
                    if self.check_toehold(e):
                        if (b.node2, b.dom2[i]) not in notbond + potbondto \
                                and(strand1, b.dom[i]) not in notbond + potbondto \
                                and self.bondgraph.check_bonded((strand1, b.dom[i])) is not None:
                            pottoehold.append((strand1, b.dom[i]))
        else:
            for b in self.bondgraph.adj[strand1]:
                if b.node2 == strand2:
                    for i in range(0, len(b.dom)):
                        e = {(strand1, b.dom[i]), (b.node2, b.dom2[i])}
                        if self.check_toehold(e):
                            if (b.node2, b.dom2[i]) not in notbond + potbondto \
                                    and (strand1, b.dom[i]) not in notbond + potbondto \
                                    and self.bondgraph.check_bonded((strand1, b.dom[i])) is not None:
                                pottoehold.append((strand1, b.dom[i]))

        pottoehold = util.get_closest_target(potbondto, pottoehold)
        if pottoehold is None:
            return None, None
        connect = self.bondgraph.get_connection(pottoehold, strand2)

        if strand1 == strand2:
            return pottoehold, pottoehold
        elif connect is None:
            return None, None
        elif len(connect) == 0:
            return None, None

        return pottoehold, connect[0]

    def have_anchor(self, startstrand, endstrand, domlow, domhigh):
        d1 = domlow - 1
        d2 = domhigh + 1
        for b in self.bondgraph.adj[startstrand]:
            if (d1 in b.dom) or (d2 in b.dom):
                if b.node2 == endstrand:
                    return -1
                else:
                    return b.node2
        return -2
    '''
        def check_adjacent_bond(self, v, n):
        """
        check if there are bonds existing adjacent to the edge
        :param v:
        :param n:
        :return:
        """
        d1 = self.bondgraph.check_bonded((v[0], n[0] - 1))
        d2 = self.bondgraph.check_bonded((v[1], n[1] - 1))
        d3 = self.bondgraph.check_bonded((v[1], n[1] + 1))
        d4 = self.bondgraph.check_bonded((v[0], n[0] + 1))

        # TODO: use another function, check correctness

        if d1 is not None:
            if d2 is not None:
                if self.bondgraph.search_path(v[0], [], v[1], self.bondgraph.get_edges()):
                    potloop = self.bondgraph.loop.pop()
                    lastdom = d1[1]
                    for i in range(1, len(potloop) - 1):
                        newbond1 = self.bondgraph.check_bonded((potloop[i], lastdom + 1))
                        newbond2 = self.bondgraph.check_bonded((potloop[i], lastdom - 1))
                        if newbond1 is not None:
                            lastdom = newbond1[1]
                            continue
                        if newbond2 is not None:
                            lastdom = newbond2[1]
                            continue
                        else:
                            return True
                    return False
                else:
                    return False

        if d4 is not None:
            if d3 is not None:
                if self.bondgraph.search_path(v[0], [], v[1], self.bondgraph.get_edges()):
                    potloop = self.bondgraph.loop.pop()
                    lastdom = d4[1]
                    for i in range(1, len(potloop) - 1):
                        newbond1 = self.bondgraph.check_bonded((potloop[i], lastdom + 1))
                        newbond2 = self.bondgraph.check_bonded((potloop[i], lastdom - 1))
                        if newbond1 is not None:
                            lastdom = newbond1[1]
                            continue
                        if newbond2 is not None:
                            lastdom = newbond2[1]
                            continue
                        else:
                            return True
                    return False
                else:
                    return False

        return True
    '''
