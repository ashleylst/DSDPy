from src.species import species_explore as se
from src.util import util, cexception
from src.basics import output as on, generate_pysbmodel as gp, initialize_system


def start_processor(filedir='../res/input', threshold=10, window=None):
    """
    the entry point to DSDPy

    :param window:
    :param threshold:
    :param filedir: file directory to the input file
    """
    # initialization
    specieslist, speciesidmap, kinetics, initnames, concentrations, outdir, simupara = \
        initialize_system.initialize(filedir)

    if outdir == '':
        outdir = '../output'
    if len(simupara) == 0:
        simupara = [100, 100]
    initlen = len(specieslist)
    if len(initnames) < initlen:
        for i in range(len(initnames), initlen):
            initnames.append('ss_'+str(i+1))
    elif len(initnames) > initlen:
        raise cexception.SpeciesError("there are more initial species names than the number of initial species")

    reactionlist = []
    visited = [False for _ in range(0, len(specieslist))]
    indexlist = []
    cursor = 0
    iteration = 0

    # explore all possibilities in species with regards to the initial DSD system
    while not visited[cursor]:
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

        if iteration == threshold:
            break
        iteration += 1

    # example use for a possible debugging option defobs :
    # md = gp.generate_model(specieslist, reactionlist, initlen, initnames, concentrations, defobs=[8, 10])

    md = gp.generate_model(specieslist, reactionlist, initlen, initnames, concentrations)

    # if there can be reactions, then simulate
    if len(md.rules) != 0:
        # example use for using Scipy ODE simulator:
        # on.simulate_scipy(md, filedir=outdir, time=simupara[0], steps=simupara[1])
        x, y, obs = on.simulate_bng(md,
                        time=simupara[0],
                        steps=simupara[1])
        on.visualize_simulation_results(x, y, obs)

    # output for GUI interface
    on.output_network_txt(specieslist,
                          reactionlist,
                          filedir=outdir)


start_processor(filedir='../res/input')


