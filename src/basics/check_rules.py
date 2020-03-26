def check_binding(graph):
    """
    check rule binding

    :param graph: a StrandGraph object
    :return: edges to be added
    """
    changeedge = []
    for i in graph.A:
        if i not in graph.E:
            if not graph.same_species(i):
                if not graph.check_toehold(i):
                    continue
            if graph.available(i):
                print("RB: " + str(i.copy()))
                changeedge.append(i)
    return changeedge


def check_unbinding(graph):
    """
    check rule unbinding

    :param graph: a StrandGraph object
    :return: edges to be deleted
    """
    changeedge = []
    for i in graph.E:
        if graph.check_toehold(i):
            if not graph.anchored(i):
                print("RU: " + str(i.copy()))
                changeedge.append(i)
    return changeedge


def check_existence_4way(elist, edge):
    """
    prevent double check in four-way migration

    :param elist: changeedge4
    :param edge: edge to be examined
    :return:
    """
    for x in elist:
        if edge == x[1]:
            return True
    return False


def check_migration(graph):
    """
    check rule migration

    :param graph: a StrandGraph object
    :return: changeedge3: edges to be changed in three-way migration
             changeedge4: edges to be changed in four-way migration
    """
    changeedge3 = []
    changeedge4 = []
    for i in graph.A:
        if i not in graph.E:
            if not check_existence_4way(changeedge4, i):
                notbond, tbr, r = graph.check_candidate(i)
                if len(tbr) != 0 and len(r) != 0:
                    if notbond is None:
                        if graph.check_switchable_2(i, tbr):
                            print("R4: " + str(i.copy()))
                            changeedge4.append((i, set(r), tbr[0], tbr[1]))
                    else:
                        if graph.check_switchable(notbond, tbr[0], i):
                            print("R3: " + str(i.copy()))
                            changeedge3.append((i, tbr[0]))
    return changeedge3, changeedge4

