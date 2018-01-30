from sympy import lcm, LM, LT, poly, expand, reduced
from sympy.abc import x, y

''' www.scholarpedia.org/article/Buchberger's_algorithm
    ----------------------------------------------------

input set of polynomials F '''

# www.lpthe.jussieu.fr/~talon/buchberger.html
setF = [x**2+2*x*y**2, x*y+2*y**3-1]

setG = setF  # initialise the Grobner basis of the ideal I
setM = []

''' set of unique pairs M from an initial set G

  G  {1,2,3,4}
  M  {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}}

  this is M = G x G excluding {Gi,Gi} and {Gi,Gj}={Gj,Gi}
'''

for i, s in enumerate(setG[:-1]):
    for t in setG[i+1:]:
        setM.append((s, t))

print(setM)

# setup complete F -> G, M

# run Buchberger's algorithm terminating when M is empty

while setM:

    # remove first pair from M
    p, q = setM[0]
    setM = setM[1:]
    print('p={}'.format(p))
    print('q={}'.format(q))
    print('M={}'.format(setM))

    # syzygy/s- polynomial
    s = expand(lcm(LM(p), LM(q))*(1/LT(p)*p - 1/LT(q)*q))

    print('s={}'.format(s))

    # divide s by all of G - otherwise known as normal form
    _, h = reduced(s, setG)

    print('h={}'.format(h))

    if h != 0:

        ''' pair h with all members of G (Cartesian/Outer product)

            P = G x h

            {1,2,3}x{4} = {{1,4},{2,4},{3,4}} '''

        # add new pairs P straight into M and h into G
        for g in setG:
            setM.append([g, h])

        print('M+P={}'.format(setM))
        setG.append(h)

print("F = {}\nG = {}".format(setF, setG))
