from sympy import lcm, LM, LT, poly, expand, reduced
from sympy.abc import x, y

def s_polynomial(f, g):
    return expand(lcm(LM(f), LM(g))*(1/LT(f)*f - 1/LT(g)*g))

def buchberger(F):
    """Toy implementation of Buchberger algorithm. """
    G, pairs = list(F), set([])

    for i, f1 in enumerate(F):
        for f2 in F[i+1:]:
            pairs.add((f1, f2))

    while pairs:
        #f1, f2 = pairs.popitem()
        print('all pairs {}'.format(pairs))
        f1, f2 = list(pairs)[0]
        pairs = set(list(pairs)[1:])

        print('trunc pairs {}'.format(pairs))

        s = s_polynomial(f1, f2)
        print('f1 {}'.format(f1))
        print('f2 {}'.format(f2))
        print('s {}'.format(s))
        _, h = reduced(s, G)
        print('h {}'.format(h))

        if h != 0:
            for g in G:
                pairs.add((g, h))

            G.append(g)

    for i, g in enumerate(G):
        _, G[i] = reduced(g, G[:i] + G[i+1:])

    G = map(monic, G)

    return G

buchberger([x*2+2*x*y**2, x*y+2*y**3-1])
