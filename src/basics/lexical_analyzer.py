from src.strand import strand as sta
from src.util import cexception as ex


def lexer_strand(str, cnt):
    """
    a lexer for transforming strand in any canonical form to a Strand object

    :param str: string
    :param cnt: count
    :return: a Strand object
    """
    state = 0

    if len(str) < 3:
        raise ex.SpeciesError("illegal species text representation.")

    strand = sta.Strand()
    for i in range(0, len(str)):
        ch = str[i]

        if i == 0:
            if ch != '<':
                raise ex.SpeciesError("illegal species text representation.")
            name = ''
            toehold = False
            comp = False
            bond = False
            bondname = ''
            continue

        if 'a' <= ch <= 'z' or 'A' <= ch <= 'Z' or '0' <= ch <= '9' and state != 5:
            if state == 0:
                name += ch
                continue
            if state == 3:
                bondname += ch
                continue

        if ch == '^':
            state = 1
            toehold = True
            continue

        if ch == '*':
            state = 2
            comp = True
            continue

        if ch == '!':
            state = 3
            bond = True
            continue

        if ch == ' ' and state != 4:
            dom = sta.Domain(name, toehold, comp, bond, bondname)
            strand.domains.append(dom)

            state = 0

            name = ''
            toehold = False
            comp = False
            bond = False
            bondname = ''
            continue

        if ch == '>':
            dom = sta.Domain(name, toehold, comp, bond, bondname)
            strand.domains.append(dom)

            state = 4

        if ch == '\n':
            if str[i-1] != '>':
                raise ex.SpeciesError("illegal species text representation.")
            strand.add_color(cnt)
            state = 0

    return strand


def lexer_site(string):
    sites = []
    name = ''

    for i in range(1, len(string) - 1):
        ch = string[i]

        if 'a' <= ch <= 'z' or 'A' <= ch <= 'Z' or '0' <= ch <= '9':
            name += ch
            continue

        if ch == '*':
            name += '_p'
            continue

        if ch == '!':
            name += '_'
            continue

        if ch == ' ':
            sites.append(name)
            name = ''
            continue

    sites.append(name)

    return sites
