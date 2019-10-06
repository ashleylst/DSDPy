class Bond:
    """
    Data structure of bonds between two domains
    """

    node1 = 0
    '''start node'''

    node2 = 0
    '''end node'''

    dom = []
    '''domain of the startnode'''

    dom2 = []
    '''domain of the endnode'''

    def __init__(self, node1, node2, dom, dom2):
        self.node1 = node1
        self.node2 = node2
        self.dom = dom
        self.dom2 = dom2

    def appenddom(self, dom, dom2):
        """
        append domain numbers

        :param dom: domain number on the start strand
        :param dom2: domain number on the end strand
        """
        self.dom.append(dom)
        self.dom2.append(dom2)

    def transform_E(self):
        """
        flatten the edges

        :return: edges in list
        """
        E = []
        for i in range(0, len(self.dom)):
            E.append(((self.node1, self.dom[i]), (self.node2, self.dom2[i])))
        return E
