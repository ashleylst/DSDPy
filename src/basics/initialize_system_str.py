from src.basics import lexical_analyzer as lex
from src.strand import strand_graph as sg, bond_graph as bg
from src.species import species as sp
from src.species import species_explore as se
from src.util import cexception as ex
from bidict import bidict


def get_additional_info(textls, pos):
    names = []
    concentrations = []
    flag = False

    for i in range(pos+1, len(textls)):
        line = textls[i]
        if line == '--':
            flag = True
            kinetics, outdir, simupara = get_kinetics(textls, i)
            break

        line = line.split(' ')
        names.append(line[0])
        concentrations.append(int(line[1]))

    if not flag:
        raise ex.KineticsError("kinetics not defined.")

    return names, concentrations, kinetics, outdir, simupara


def get_kinetics(textls, pos):

    kinetics = {}
    outdir = ''
    simupara = []
    for i in range(pos+1, len(textls)):
        line = textls[i]
        if line == '--':
            outdir, simupara = get_outdir_simupara(textls, i)
            break

        line = line.split(' ')
        kinetics[line[0]] = float(line[1])
    return kinetics, outdir, simupara


def get_outdir_simupara(textls, pos):

    outdir = ''
    simupara = [100, 100]
    for i in range(pos+1, len(textls)):
        line = textls[i]
        line = line.split('=')
        if line[0] == 'outdir':
            outdir = line[1]
            continue
        if line[0] == 'time':
            simupara[0] = line[1]
        if line[0] == 'time_step':
            simupara[1] = line[1]
    return outdir, simupara


def initialize(text):

    kinetics = {}
    strands = []
    speciesbreak = []
    splabelling = []
    speciesnum = 1
    cnt = 0
    cursplabel = ''
    textls = text.split('\n')
    while '' in textls:
        textls.remove('')

    for i in range(len(textls)):
        line = textls[i]

        if cursplabel == '':
            cursplabel += line[0:len(line)]
        else:
            cursplabel += '|' + line[0:len(line)]

        if line == '--':
            splabelling.append(cursplabel[0:len(cursplabel)-3])
            names, concentrations, kinetics, outdir, simupara = get_additional_info(textls, i)
            break

        if line == '//':
            speciesbreak.append(len(strands))
            splabelling.append(cursplabel[0:len(cursplabel)-3])
            cursplabel = ''
            speciesnum += 1
            continue

        strand = lex.lexer_strand(line, cnt)
        strand.add_color(cnt)
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
