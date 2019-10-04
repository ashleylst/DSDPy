class Strand:
    """
    Data structure of the strand
    """

    domains = []
    '''domains of the strand'''

    color = 0
    '''color of the strand'''

    def __init__(self):
        self.domains = []

    def add_domain(self, domain):
        self.domains.append(domain)

    def add_color(self, color):
        self.color = color

    def check_same_strand(self, strand1):
        """
        check if strand1 is equivalent to the current strand

        :param strand1:
        :return: True if equivalent, False otherwise
        """
        if len(self.domains) != len(strand1.domains):
            return False
        else:
            for i in range(0, len(self.domains)):
                if self.domains[i].check_same_domain(strand1.domains[i]):
                    continue
                else:
                    return False
            return True


class Domain:
    """
    Data structure of the domain
    """

    name = ''
    '''name of the domain'''

    toehold = False
    '''boolean variable indicating if the domain is a toehold'''

    comp = False
    '''boolean variable indicating if the domain is a complementary domain to some other domain'''

    bond = False
    '''boolean variable indicating if the domain is bonded'''

    bondname = ''
    '''name of the bond, empty if not bonded'''

    def __init__(self, name, toehold, comp, bond, bondname):
        self.name = name
        self.toehold = toehold
        self.comp = comp
        self.bond = bond
        self.bondname = bondname

    def set_bond(self, bondname):
        self.bond = True
        self.bondname = bondname

    def check_same_domain(self, domain1):
        """
        check if domain1 is equivalent to the current domain

        :param domain1:
        :return: True if equivalent, False otherwise
        """
        if self.name == domain1.name and self.toehold == domain1.toehold and self.comp == domain1.comp:
            return True
        return False
