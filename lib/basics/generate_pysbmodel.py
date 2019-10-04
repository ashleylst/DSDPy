from pysb import *


def generate_monomer(species, monomerdict, initlen, initnames, tbobs):
    """
    generate a PySB monomer based on species

    :param species: a Species object
    :param monomerdict: a dictionary with all monomers linked to their species id
    :param initlen: number of the initial species
    :param initnames: names of the initial species
    :return: monomer, observable, and updated monomerdict
    """
    if species.id <= initlen:
        name = initnames[species.id - 1]
    else:
        name = 'sp_' + str(species.id)
    sites = ['init']

    m = Monomer(name, sites)
    monomerdict[species.id] = m
    if species.id in tbobs:
        Observable('obs' + name, m(init=None))
    return monomerdict


def generate_rule(reaction, i, monomerdict, rate):
    """
    generate a PySB rule based on reaction

    :param reaction: a Reaction object
    :param i: ordered position of reaction in reactionlist
    :param monomerdict: a dictionary with all monomers linked to their species id
    :param rate: reaction rate, a Parameter object
    :return: PySB rule
    """
    name = reaction.rule + '_' + str(i)

    reactants = []
    for i in reaction.reactants:
        r = monomerdict[i.id]
        reactants.append(r)

    products = []
    for i in reaction.products:
        p = monomerdict[i.id]
        products.append(p)

    if len(reactants) == 1:
        if len(products) == 1:
            rule = Rule(name, reactants[0](init=None) >>
                        products[0](init=None), rate)
        else:
            rule = Rule(name, reactants[0](init=None) >>
                        products[0](init=None) + products[1](init=None), rate)
    else:
        if len(products) == 1:
            rule = Rule(name, reactants[0](init=None) + reactants[1](init=None) >>
                        products[0](init=None), rate)
        else:
            rule = Rule(name, reactants[0](init=None) + reactants[1](init=None) >>
                        products[0](init=None) + products[1](init=None), rate)

    return rule


def generate_concentrations(initlen, concentrations):
    """
    generate concentrations of the initial species

    :param initlen: number of the initial species
    :param concentrations: concentrations of the initial species
    :return: PySB parameters for concentrations
    """
    para = []
    for i in range(0, initlen):
        name = 'ic' + str(i)
        para.append(Parameter(name, concentrations[i]))
    return para


def generate_model(specieslist, reactionlist, initlen, initnames, concentrations, defobs=list()):
    """
    generate a PySB model for the DSD reaction network

    :param specieslist: list of species that can be synthesized in the DSD system
    :param reactionlist: list of reactions that can happen in the DSD system
    :param initlen: number of the initial species
    :param initnames: names of the initial species
    :param concentrations: concentrations of the initial species
    :param defobs: a debug argument for defining observables
    :return: a PySB model
    """
    model = Model()

    monomerdict = {}

    if len(defobs) == 0:
        tbobs = [i + 1 for i in range(0, len(specieslist))]
    else:
        tbobs = [i + 1 for i in range(0, initlen)] + defobs

    for i in specieslist:
        monomerdict = generate_monomer(i, monomerdict, initlen, initnames, tbobs)
        # model.monomers.add(monomer)
        # model.observables.add(obs)

    count = 0
    for i in reactionlist:
        count += 1
        rate = Parameter('k_' + str(count), i.rate)
        model.rules.add(generate_rule(i, count, monomerdict, rate))

    # Parameter('ic', 10000)
    para = generate_concentrations(initlen, concentrations)
    for i in range(0, initlen):
        Initial(monomerdict[i+1](init=None), para[i])

    return model
