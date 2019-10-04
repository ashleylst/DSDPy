from lib.basics import lexical_analyzer as lex
from lib.strand import strand_graph as sg, bond_graph as bg
from lib.species import species as sp
from lib.species import species_explore as se
from bidict import bidict


def get_additional_info(fp, line):
    """
    get additional info from input

    :param fp: opened file
    :param line: current line
    :return: names of initial species,
            concentrations of initial species,
            a dictionary object of kinetics,
            output directory,
            simulation parameters
    """
    names = []
    concentrations = []
    while line:
        line = fp.readline()
        if line == '--\n':
            kinetics, outdir, simupara = get_kinetics(fp, line)
            break
        line = line.strip('\n')
        line = line.split(' ')
        names.append(line[0])
        concentrations.append(int(line[1]))
    return names, concentrations, kinetics, outdir, simupara


def get_kinetics(fp, line):
    """
    get the kinetics info from input

    :param fp: opened file
    :param line: current line
    :return: a dictionary object of kinetics,
            output directory,
            simulation parameters
    """
    kinetics = {}
    outdir = ''
    simupara = []
    while line:
        line = fp.readline()
        if not line:
            break
        if line == '--\n':
            outdir, simupara = get_outdir_simupara(fp, line)
            break
        line = line.strip('\n')
        line = line.split(' ')
        kinetics[line[0]] = float(line[1])
    return kinetics, outdir, simupara


def get_outdir_simupara(fp, line):
    """
    get output directory and simulation parameters (in format: time steps)

    :param fp: opened file
    :param line: current line
    :return: output directory and simulation  parameters
    Note: output directory is '' if not specified and
    list of simulation parameters is empty if not specified
    """
    outdir = ''
    simupara = []
    while line:
        line = fp.readline()
        if not line:
            break
        if line == '--\n':
            continue
        line = line.strip('\n')
        line = line.split(' ')
        if len(line) == 1:
            outdir = line
        else:
            simupara = line
    return outdir, simupara


def initialize(filedir):
    """
    initialize the DSD system

    :param filedir: file directory of the input
    :return: specieslist,
            speciesidmap,
            kinetics,
            names,
            concentrations,
            output directory,
            simulation parameters
    """
    kinetics = {}
    strands = []
    speciesbreak = []

    with open(filedir) as fp:
        line = fp.readline()
        cnt = 1

        strand = lex.lexer_strand(line, cnt)
        strands.append(strand)
        cnt += 1

        while line:
            # print("Line {}: {}".format(cnt, line.strip()))
            line = fp.readline()

            if not line:
                break

            if line == '--\n':
                names, concentrations, kinetics, outdir, simupara = get_additional_info(fp, line)
                break

            if line == '//\n':
                speciesbreak.append(len(strands))
                continue

            strand = lex.lexer_strand(line, cnt)
            for i in range(0, len(strands)):
                if strands[i].check_same_strand(strand):
                    strand.add_color(strands[i].color)
                    cnt -= 1
                    break
            strands.append(strand)
            cnt += 1

    specieslist = []
    speciesidmap = bidict()

    for i in range(0, len(speciesbreak) + 1):
        if i == 0:
            low = 0
        else:
            low = speciesbreak[i - 1]
        if i == len(speciesbreak):
            high = len(strands)
        else:
            high = speciesbreak[i]
        strandgraph = sg.StrandGraph(strands[low:high])
        colorset, colormap = se.generate_colorinfo(strandgraph.color)

        species = sp.Species(strandgraph.V, colorset, colormap, strandgraph)
        species.set_id(i + 1)
        speciesidmap.put(species.id, species.canonicalform)
        specieslist.append(species)

    '''
    strandgraph = sg.StrandGraph(strands)

    speciesnodes = strandgraph.bondgraph.get_species()
    specieslist = []

    cnt = 0
    speciesidmap = bidict()

    for i in speciesnodes:
        cnt += 1

        sub = bg.SubBondGraph(i, strandgraph.color, strandgraph.bondgraph.adj, strandgraph.V)

        species = sp.Species(i, sub.colorset, sub.colormap, strandgraph)
        species.set_id(cnt)

        speciesidmap.put(species.id, species.canonicalform)

        specieslist.append(species)
    '''
    return specieslist, speciesidmap, kinetics, names, concentrations, outdir, simupara
