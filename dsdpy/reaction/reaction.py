class Reaction:
    """
    Data structure for a reaction in DSD system
    """
    reactants = []

    products = []

    rule = ''
    '''the rule this reaction applies'''

    rate = -1
    '''reaction rate'''

    def __init__(self, reactants, products):
        self.reactants = reactants
        self.products = products

    def add_rate(self, rate):
        self.rate = rate

    def add_product(self, product):
        self.products.append(product)

    def add_rule(self, rulename):
        self.rule = rulename

    def generate_output(self):
        """
        generate string output for output txt file
        :return: string representation of the reaction. e.g. 1 + 2 --> 3 rate=1
        """
        string = self.rule + ' '
        rlen = len(self.reactants)
        plen = len(self.products)

        for i in range(0, rlen):
            if i != rlen - 1 and rlen != 1:
                string += str(self.reactants[i].id) + ' + '
            else:
                string += str(self.reactants[i].id)

        string += ' --> '

        for i in range(0, plen):
            if i != plen - 1 and plen != 1:
                string += str(self.products[i].id) + ' + '
            else:
                string += str(self.products[i].id)

        string += ' rate=' + str(self.rate) + '\n'

        return string
