from src.strand import strand_graph as sg, bond_graph as bg
from src.basics import check_rules as cr
from src.species import species as sp
from src.util import util
from src.reaction import reaction as ra
from collections import defaultdict
import copy


def check_existence(species, specieslist, speciesidmap):
    """
    check if the produced species is new in the DSD system

    :param species: produced species
    :param specieslist: list of species in the system
    :param speciesidmap: bi-directional map of species id and species canonical form
    :return: specieslist, speciesidmap
    """
    if not speciesidmap.inverse.__contains__(species.canonicalform):
        species.set_id(len(specieslist) + 1)
        specieslist.append(species)
        speciesidmap.put(species.id, species.canonicalform)
        pos = species.id
    else:
        pos = speciesidmap.inverse[species.canonicalform]
    return specieslist, speciesidmap, pos


def generate_colorinfo(color):
    """
    generate color information based on strand graph's color

    :param color: strand graph's color
    :return: colorset: a set of color
             colormap: a map of color to id of the strand
    """
    colorset = set()
    colormap = defaultdict(set)

    for i in range(0, len(color)):
        colorset.add(color[i])
        colormap[color[i]].add(i)
    return colorset, colormap


def generate_species(strandgraph, specieslist, speciesidmap, reaction):
    """
    generate a species

    :param reaction:
    :param strandgraph:
    :param specieslist:
    :param speciesidmap:
    :return: specieslist, speciesidmap
    """
    colorset, colormap = generate_colorinfo(strandgraph.color)
    newspecies = sp.Species(strandgraph.V, colorset, colormap, strandgraph)
    specieslist, speciesidmap, pos = check_existence(newspecies, specieslist, speciesidmap)
    reaction.add_product(specieslist[pos - 1])
    return specieslist, speciesidmap, reaction


def generate_multiple_species(strandgraph, specieslist, speciesidmap, reaction):
    """
    generate multiple species (applies when unbinding and migrating)

    :param reaction:
    :param strandgraph:
    :param specieslist:
    :param speciesidmap:
    :return:
    """
    if strandgraph.bondgraph.speciesnum != 1:
        speciesnodes = strandgraph.bondgraph.get_species()
        for i in speciesnodes:
            sub = bg.SubBondGraph(i, strandgraph.color, strandgraph.bondgraph.adj, strandgraph.V)
            newspecies = sp.Species(i, sub.colorset, sub.colormap, strandgraph)
            specieslist, speciesidmap, pos = check_existence(newspecies, specieslist, speciesidmap)
            reaction.add_product(specieslist[pos - 1])
    else:
        specieslist, speciesidmap, reaction = generate_species(strandgraph, specieslist, speciesidmap, reaction)
    return specieslist, speciesidmap, reaction


def check_hidden_prevbond(E, strandgraph):
    hiddenset = set(strandgraph.bondgraph.hidden)
    for e in E:
        if hiddenset & e:
            return True
    return False


def mono(species, specieslist, speciesidmap, reactionlist, kinetics):
    """
    mono-molecule reaction

    :param reactionlist:
    :param species:
    :param specieslist:
    :param speciesidmap:
    :return: specieslist, speciesidmap
    """
    strandgraph = sg.StrandGraph(species.strands)

    prevE = strandgraph.E
    prevSG = copy.copy(strandgraph)
    flag = False

    e3, e4 = cr.check_migration(prevSG)
    if len(e3) != 0:
        miggroup = strandgraph.bondgraph.check_following_migration(e3)

        #miggroup = util.check_following_migration(e3, p=0)

        for x in miggroup:
            startv = list(e3[x[0]][0] - e3[x[0]][1])
            sortdom = util.get_migrate_nodes(e3, x, startv[0][0])
            '''
            endv = list(e3[x[0]][0] & e3[x[0]][1])
            p = prevSG.have_anchor(startv[0][0], endv[0][0], sortdom[0], sortdom[len(sortdom)-1])
            if p != -1:
                if p == -2:
                    continue
                elif p >= 0:
                    con = prevSG.bondgraph.get_connection((p, 0), endv[0][0])
                    if len(con) == 0:
                        continue
            '''
            reaction = ra.Reaction([species], [])
            E = copy.copy(prevE)
            notbond = []
            potbondto = []
            for i in range(0, len(x)):
                if not strandgraph.check_hidden_all(e3[x[i]]):
                    flag = True
                E.remove(e3[x[i]][1])
                E.append(e3[x[i]][0])
                notbond += list(e3[x[i]][0] - e3[x[i]][1])
                potbondto += list(e3[x[i]][0] & e3[x[i]][1])

            '''
            target1, target2 = strandgraph.get_connect_toehold(notbond, potbondto)
            if target1 is None or target2 is None:
                    continue

            closest1, closest2 = util.get_domains_on_2sides(target1, target2, notbond, potbondto)
            if not strandgraph.check_complementary(closest1, closest2):
                continue
            '''
            if not flag:
                continue

            strandgraph.reconstruct(E)

            anchor = False
            for j in range(len(x)):
                if (startv[0][0], sortdom[0]) in e3[x[j]][0] or (startv[0][0], sortdom[len(sortdom)-1]) in e3[x[j]][0]:
                    if strandgraph.anchored(e3[x[j]][0]):
                        anchor = True
                        break
            if not anchor:
                strandgraph.reconstruct(prevE)
                continue

            specieslist, speciesidmap, reaction = \
                generate_multiple_species(strandgraph, specieslist, speciesidmap, reaction)
            reaction.add_rule('R3')
            reaction.add_rate(kinetics['R3'])
            reactionlist.append(reaction)

    if len(e4) != 0:
        for x in e4:
            reaction1 = ra.Reaction([species], [])
            E = copy.copy(prevE)
            E.remove(x[2])
            E.remove(x[3])
            E.append(x[0])
            E.append(x[1])
            strandgraph.reconstruct(E)

            specieslist, speciesidmap, reaction1 = \
                generate_multiple_species(strandgraph, specieslist, speciesidmap, reaction1)
            reaction1.add_rule('R4')
            reaction1.add_rate(kinetics['R4'])
            reactionlist.append(reaction1)

    e = cr.check_binding(prevSG)
    if len(e) != 0:
        for x in e:
            reaction2 = ra.Reaction([species], [])
            E = copy.copy(prevE)
            E.append(x)
            strandgraph.reconstruct(E)

            # check if there is a previous bond hidden
            # if check_hidden_prevbond(prevE, strandgraph):
            #   continue

            specieslist, speciesidmap, reaction2 = \
                generate_species(strandgraph, specieslist, speciesidmap, reaction2)
            reaction2.add_rule('RB')
            reaction2.add_rate(kinetics['RB'])
            reactionlist.append(reaction2)

    e = cr.check_unbinding(prevSG)
    if len(e) != 0:
        for x in e:
            reaction3 = ra.Reaction([species], [])
            E = copy.copy(prevE)
            E.remove(x)
            strandgraph.reconstruct(E)

            specieslist, speciesidmap, reaction3 = \
                generate_multiple_species(strandgraph, specieslist, speciesidmap, reaction3)
            reaction3.add_rule('RU')
            reaction3.add_rate(kinetics['RU'])
            reactionlist.append(reaction3)

    return specieslist, speciesidmap, reactionlist


def bi(speciescomb, specieslist, speciesidmap, reactionlist, kinetics):
    """
    bi-molecule reaction

    :param reactionlist:
    :param speciescomb: combinations of species ids
    :param specieslist:
    :param speciesidmap:
    :return: specieslist, speciesidmap
    """
    index1, index2 = speciescomb
    species1 = specieslist[index1]
    species2 = specieslist[index2]
    len1 = len(species1.nodes)

    strandgraph = sg.StrandGraph(species1.strands + species2.strands, merge=len1)

    prevE = strandgraph.E
    prevSG = copy.copy(strandgraph)

    e = cr.check_binding(prevSG)
    if len(e) != 0:
        for x in e:
            reaction = ra.Reaction([species1, species2], [])
            # check if the bond is between two species
            v, n = util.get_edge_info(x)
            if (v[0] < len1 and v[1] < len1) or (v[0] >= len1 and v[1] >= len1):
                continue

            E = copy.copy(prevE)
            E.append(x)
            strandgraph.reconstruct(E)

            # if check_hidden_prevbond(prevE, strandgraph):
            #    continue

            specieslist, speciesidmap, reaction = generate_species(strandgraph, specieslist, speciesidmap, reaction)
            reaction.add_rule('RB')
            reaction.add_rate(kinetics['RB'])
            reactionlist.append(reaction)

    return specieslist, speciesidmap, reactionlist
