import networkx as nx
import numpy as np
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
from pysb.bng import *


def generate_incidence_matrix(specieslist, reactionlist):
    """

    :param specieslist: list of species
    :param reactionlist: list of reactions
    :return:
    """
    nodes = [i.id for i in specieslist]
    edges = []

    for i in reactionlist:
        for j in i.reactants:
            for k in i.products:
                edges.append([j.id, k.id])

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    matrix = -nx.incidence_matrix(G, oriented=True)

    return nodes, edges, matrix


def visualize_simulation_results(x, y, obs, filedir, option='bng'):
    """

    :param x: values on x axis
    :param y: values on y axis
    :param obslen: length of the observables
    :param initlen: length of the initial species
    :param initnames: names of the initial species
    :param filedir: file directory to write the image
    """
    plt.figure()

    obslen = len(obs)

    for i in range(0, obslen):
        label = obs[i].name[3:]
        if option == 'bng':
            plt.plot(x, y[:, i], label=label)
        elif option == 'scipy':
            plt.plot(x, y[obs[i].name], label=label)

    plt.xlabel("Time (s)")
    plt.ylabel("Complexes")
    plt.legend(bbox_to_anchor=(1.04, 1))

    plt.savefig(filedir + '/simres.png', bbox_inches='tight', pad_inches=0.5)


def output_network_txt(specieslist, reactionlist, filedir='../output'):
    """
    Text file for GUI to visualize reaction network

    :param specieslist: list of species
    :param reactionlist: list of reactions
    :param filedir: file directory to write the txt file
    """
    file = open(filedir + '/output.txt', 'w+')
    file.write('-----Species-----\n')
    for i in specieslist:
        file.write(i.generate_output())
        file.write('\n')

    file.write('-----Reactions-----\n')
    for i in reactionlist:
        file.write(i.generate_output())
        file.write('\n')

    file.write('-----Incidence Matrix-----\n')
    rowlabels, collabels, incidencematrix = generate_incidence_matrix(specieslist, reactionlist)
    incidencematrix = incidencematrix.todense()

    print(*collabels, file=file)
    np.set_printoptions(linewidth=90)

    for rowlabel, row in zip(rowlabels, incidencematrix):
        print('%s %s' % ('%03s' % rowlabel, ' '.join('%s' % i for i in row)), file=file)

    file.close()


def simulate_BNG(model, time=1000, steps=100, bngnetwork=False, filedir='../output'):
    """
    simulate the reaction network using BNG

    :param model: a PySB object
    :param initlen: length of initial species
    :param initnames: names of initial species
    :param time: simulation time
    :param steps: simulation steps
    :param bngnetwork: a boolean variable indicating if BNG output file is needed
    :param filedir: file directory to write output files
    """
    # TODO: generate BNG output file network
    if bngnetwork:
        network = generate_network(model=model)

    obslen = len(model.observables)

    output = run_ssa(model=model, t_end=time, n_steps=steps)
    output = output.tolist()
    output = np.array(output)

    visualize_simulation_results(output[:, 0], output[:, obslen + 1: obslen * 2 + 1], model.observables, filedir, option='bng')


def simulate_Scipy(model, time=1000, steps=100, filedir='../output'):
    """
    simulate the reaction network using Scipy ODE

    :param model: a PySB object
    :param time: simulation time
    :param steps: simulation steps
    :param filedir: file directory to write output files
    """
    t = np.linspace(0, time, steps)
    simres = ScipyOdeSimulator(model, tspan=t).run()
    yout = simres.all

    visualize_simulation_results(t, yout, model.observables, filedir, option='scipy')