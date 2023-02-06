from src.basics import lexical_analyzer as lex
from src.strand import strand_graph as sg, bond_graph as bg
from src.species import species as sp
from src.species import species_explore as se
from src.util import cexception as ex
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
    flag = False

    while line:
        line = fp.readline()
        if line == '--\n':
            flag = True
            kinetics, outdir, simupara = get_kinetics(fp, line)
            break
        if not line:
            break
        line = line.strip('\n')
        line = line.split(' ')
        names.append(line[0])
        concentrations.append(int(line[1]))

    if not flag:
        raise ex.KineticsError("kinetics not defined.")

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
            outdir = str(line[0])
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
    splabelling = []
    speciesnum = 0

    with open(filedir) as fp:
        line = fp.readline()
        cnt = 1

        strand = lex.lexer_strand(line, cnt)
        strands.append(strand)
        speciesnum += 1
        cursplabel = line[0:len(line)-1]
        cnt += 1

        while line:
            # print("Line {}: {}".format(cnt, line.strip()))
            line = fp.readline()
            if cursplabel == '':
                cursplabel += line[0:len(line)-1]
            else:
                cursplabel += '|' + line[0:len(line)-1]

            if not line:
                break

            if line == '--\n':
                splabelling.append(cursplabel[0:len(cursplabel)-3])
                names, concentrations, kinetics, outdir, simupara = get_additional_info(fp, line)
                break

            if line == '//\n':
                speciesbreak.append(len(strands))
                splabelling.append(cursplabel[0:len(cursplabel)-3])
                cursplabel = ''
                speciesnum += 1
                continue

            strand = lex.lexer_strand(line, cnt)
            for i in range(0, len(strands)):
                if strands[i].check_same_strand(strand):
                    strand.add_color(strands[i].color)
                    cnt -= 1
                    break
            strands.append(strand)
            cnt += 1

    if speciesnum == len(splabelling) + 1:
        splabelling.append(cursplabel)

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

    # error handling
    try:
        if len(kinetics) < 4:
            raise ex.KineticsError("not sufficient types of reaction rates defined.")
    except NameError:
        raise ex.KineticsError("reaction rates undefined.")

    for r in ['RB', 'RU', 'R3', 'R4', 'RB_2']:
        if kinetics[r] == 0.:
            print(r + " is set to 0.")
            raise KeyError(r + " is not defined correctly")



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
