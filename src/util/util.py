from itertools import combinations
import copy


def get_reverse(n):
    if n == 1:
        return 0
    else:
        return 1


def get_edge_info(e):
    v = [0 for i in range(2)]
    n = [0 for i in range(2)]
    t = 0

    for x in e:
        v[t], n[t] = x
        t += 1
    return v, n


def sort_e_by_domain(val):
    return val[0][1]


def sort_by_strand(val):
    return val[0][0]


def check_edge_in_tuplelist(edge, tpl):
    for i in tpl:
        if edge in i:
            return True
    return False


def compare(a, b):
    return (a > b) - (a < b)


def flip(i):
    if i == 0:
        i = 1
    elif i == 1:
        i = 0
    return i


def get_free_domains(limits, blocks, bound):
    limits = sorted(limits)
    interval = limits[1] - limits[0]
    for i in blocks:
        if limits[1] > i > limits[0]:
            tmp = abs(bound - i)
            if tmp < interval:
                interval = tmp
    return interval


def get_combinations(oldlen, newlen, cursor, indexlist):
    combold = list(combinations(indexlist[cursor:oldlen], 2))
    combself = [(i, i) for i in range(0, oldlen)]
    combnew = []

    if oldlen != newlen:
        for i in range(0, oldlen):
            for j in range(oldlen, newlen):
                combnew.append((i, j))
    return combold + combnew + combself


def get_migrate_nodes(edges, indices, startstrand):
    d = []
    for i in indices:
        vi, ni = get_edge_info(edges[i][0])
        if vi[0] == startstrand:
            d.append(ni[0])
        else:
            d.append(ni[1])
    d.sort()
    return d


def check_following_migration(edges):
    """

    :param edges:
    :return:
    """
    e = copy.copy(edges)
    visited = [False for _ in e]
    miggroup = []
    cnt = -1
    for i in range(0, len(e)):
        if visited[i]:
            continue
        e[i] = list(e[i])
        e[i][0] = list(e[i][0])
        t1 = sorted(e[i][0], key=lambda tup: tup[0])

        if not visited[i]:
            visited[i] = True
            miggroup.append([i])
            cnt += 1

        for j in range(0, len(e)):
            if j != i and not visited[j]:
                e[j] = list(e[j])
                e[j][0] = list(e[j][0])
                t2 = sorted(e[j][0], key=lambda tup: tup[0])
                if (t2[0][0] != t1[0][0]) or (t2[1][0] != t1[1][0]):
                    continue
                for num in range(0, len(miggroup[cnt])):
                    t1 = sorted(e[miggroup[cnt][num]][0], key=lambda tup: tup[0])
                    if (t1[0][1] + 1 == t2[0][1] and t1[1][1] - 1 == t2[1][1]) \
                            or (t1[0][1] - 1 == t2[0][1] and t1[1][1] + 1 == t2[1][1]):
                        visited[j] = True
                        miggroup[cnt].append(j)
                        break
    return miggroup


def get_absdist(domain1, domain2):
    """

    :param domain1:
    :param domain2:
    :return:
    """
    return abs(domain1[1] - domain2[1])


def get_closet_domain_to_target(target, domains):
    """

    :param target:
    :param domains:
    :return:
    """
    closet = 10000
    closetd = ()
    for i in domains:
        dist = get_absdist(i, target)
        if dist < closet:
            closet = dist
            closetd = i
    return closetd


def get_domains_on_2sides(target1, target2, domains1, domains2):
    """

    :param target1:
    :param target2:
    :param domains1:
    :param domains2:
    :return:
    """
    if target1[0] == domains1[0][0]:
        closetd1 = get_closet_domain_to_target(target1, domains1)
    elif target2[0] == domains1[0][0]:
        closetd1 = get_closet_domain_to_target(target2, domains1)
    if target1[0] == domains2[0][0]:
        closetd2 = get_closet_domain_to_target(target1, domains2)
    elif target2[0] == domains2[0][0]:
        closetd2 = get_closet_domain_to_target(target2, domains2)
    return closetd1, closetd2


def get_closest_target(domains, targets):
    """

    :return:
    """
    domains = sorted(domains, key=lambda tup: tup[1])
    mindist = 10000
    mint = None
    for t in targets:
        dist = min(get_absdist(t, domains[0]), get_absdist(t, domains[len(domains) - 1]))
        if dist < mindist:
            mint = t
    return mint


def check_continuity(a, b):
    for i in a:
        for j in b:
            if i + 1 == j or i - 1 == j:
                return True
    return False

