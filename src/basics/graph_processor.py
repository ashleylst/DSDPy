from src.species import species_explore as se
from src.util import util, cexception
from src.basics import output as on, generate_pysbmodel as gp, initialize_system,initialize_system_str


def initiation(filedir=None, text=None):
    if filedir:
        specieslist, speciesidmap, kinetics, initnames, concentrations, outdir, simupara = \
            initialize_system.initialize(filedir)
    elif text:
        specieslist, speciesidmap, kinetics, initnames, concentrations, outdir, simupara = \
            initialize_system_str.initialize(text)
    else:
        # TODO: Define error
        return

    if outdir == '':
        outdir = '../output'
    if len(simupara) == 0:
        simupara = [100, 100]
    initlen = len(specieslist)
    if len(initnames) < initlen:
        for i in range(len(initnames), initlen):
            initnames.append('ss_' + str(i + 1))
    elif len(initnames) > initlen:
        raise cexception.SpeciesError("there are more initial species names than the number of initial species")

    reactionlist = []
    visited = [False for _ in range(0, len(specieslist))]
    indexlist = []
    cursor = 0

    return (specieslist, speciesidmap, reactionlist, kinetics, indexlist, cursor, visited), initnames, concentrations, outdir, simupara, initlen


def one_iteration(specieslist, speciesidmap, reactionlist, kinetics, indexlist, cursor, visited):
    indexlist += [i for i in range(cursor, len(specieslist))]
    oldlen = len(visited)

    for i in range(cursor, oldlen):
        specieslist, speciesidmap, reactionlist = se.mono(specieslist[i],
                                                          specieslist,
                                                          speciesidmap,
                                                          reactionlist,
                                                          kinetics)
        visited[i] = True

    newlen = len(specieslist)

    comb = util.get_combinations(oldlen, newlen, cursor, indexlist)

    for i in comb:
        specieslist, speciesidmap, reactionlist = se.bi(i,
                                                        specieslist,
                                                        speciesidmap,
                                                        reactionlist,
                                                        kinetics)

    if oldlen != len(specieslist):
        cursor = oldlen
    else:
        cursor = oldlen - 1
    for i in range(oldlen, len(specieslist)):
        visited.append(False)

    return specieslist, speciesidmap, reactionlist, kinetics, indexlist, cursor, visited


def post_enumeration(specieslist, reactionlist, initlen, initnames, concentrations, outdir, simupara):
    md = gp.generate_model(specieslist, reactionlist, initlen, initnames, concentrations)

    if len(md.rules) != 0:
        # example use for using Scipy ODE simulator:
        # on.simulate_scipy(md, filedir=outdir, time=simupara[0], steps=simupara[1])
        x, y, obs = on.simulate_bng(md, time=simupara[0], steps=simupara[1])
    else:
        return

    text = on.generate_text(specieslist, reactionlist)

    return x, y, obs, text


def entry(filedir):
    """
    For debugging use.
    :param filedir: input file directory
    :return:
    """
    info, initnames, concentrations, outdir, simupara, initlen = initiation(filedir)
    while not info[6][info[5]]:
        info = one_iteration(*info)
    post_enumeration(info[0], info[2], initlen, initnames, concentrations, outdir, simupara)



