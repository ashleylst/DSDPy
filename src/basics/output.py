import networkx as nx
import numpy as np
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
from pysb.bng import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg


class Canvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = plt.Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.axes = fig.add_subplot(111, xmargin=0.01)
        super(Canvas, self).__init__(fig)


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


def visualize_simulation_results(x, y, obs, filedir='../output', option='bng', colormap='Paired'):
    """

    :param colormap: colormap setting
    :param option: choose visualization for bng or scipy simulation
    :param obs: observables
    :param x: values on x axis
    :param y: values on y axis
    :param filedir: file directory to write the image
    """
    plt.figure()

    obslen = len(obs)

    cmap = plt.cm.get_cmap(colormap, obslen)

    for i in range(0, obslen):
        label = obs[i].name[3:]
        if option == 'bng':
            plt.plot(x, y[:, i], label=label, c=cmap(i))
        elif option == 'scipy':
            plt.plot(x, y[obs[i].name], label=label, c=cmap(i))

    plt.xlabel("Time (s)")
    plt.ylabel("Complexes")
    plt.legend(bbox_to_anchor=(1.04, 1))

    plt.savefig(filedir + '/simres', bbox_inches='tight', pad_inches=0.5)


def output_network_txt(specieslist, reactionlist, filedir='../output'):
    """
    Text file for GUI to visualize reaction network

    :param specieslist: list of species
    :param reactionlist: list of reactions
    :param filedir: file directory to write the txt file
    """
    file = open(filedir + '/output.txt', 'w+')
    '''
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
    '''
    file.write(generate_text(specieslist, reactionlist))
    file.close()


def generate_text(specieslist, reactionlist):
    text = '-----Species-----\n'
    for i in specieslist:
        text += i.generate_output() + '\n'

    text += '-----Reactions-----\n'
    for i in reactionlist:
        text += i.generate_output() + '\n'

    text += '-----Incidence Matrix-----\n'
    rowlabels, collabels, incidencematrix = generate_incidence_matrix(specieslist, reactionlist)
    incidencematrix = incidencematrix.todense()

    text += ' '.join([str(elem) for elem in collabels]) + '\n'
    for rowlabel, row in zip(rowlabels, incidencematrix):
        text += '%s %s' % ('%03s' % rowlabel, ' '.join('%s' % i for i in row)) + '\n'

    return text


def simulate_bng(model, time=1000, steps=100, bngnetwork=False):
    """
    simulate the reaction network using BNG

    :param colormap: colormap setting
    :param model: a PySB object
    :param time: simulation time
    :param steps: simulation steps
    :param bngnetwork: a boolean variable indicating if BNG output file is needed
    :param filedir: file directory to write output files
    """
    # TODO: generate BNG output file network
    if bngnetwork:
        network = generate_network(model=model)

    monomerlen = len(model.monomers)

    output = run_ssa(model=model, t_end=time, n_steps=steps)
    output = output.tolist()
    output = np.array(output)
    row, column = output.shape

    '''
    visualize_simulation_results(output[:, 0],
                                 output[:, monomerlen + 1: column],
                                 model.observables, filedir,
                                 option='bng',
                                 colormap=colormap)
    '''

    return output[:, 0], output[:, monomerlen + 1: column], model.observables


def simulate_scipy(model, time=1000, steps=100, filedir='../output', colormap='Paired'):
    """
    simulate the reaction network using Scipy ODE

    :param colormap: colormap setting
    :param model: a PySB object
    :param time: simulation time
    :param steps: simulation steps
    :param filedir: file directory to write output files
    """
    t = np.linspace(0, time, steps)
    simres = ScipyOdeSimulator(model, tspan=t).run()
    yout = simres.all

    visualize_simulation_results(t,
                                 yout,
                                 model.observables,
                                 filedir,
                                 option='scipy',
                                 colormap=colormap)
