from math import gcd

def dotp(a, b): return sum([x * y for x, y in zip(a,b)], start=Fraction(0))

class Fraction():
    def __init__(self, num, den=1):
        self.num = num
        self.den = den
    
    def reduce(self):
        g = gcd(self.num, self.den)
        if g != 1:
            #print(g, f'{self.num}/{self.den} reduced to {self.num // g}/{self.den // g}')
            self.num //= g
            self.den //= g
        return self

    def __repr__(self): return f"({self.num}/{self.den})"

    def __add__(self, f):
        return Fraction(self.num * f.den + self.den * f.num, self.den * f.den).reduce()

    def __sub__(self, f):
        return Fraction(self.num * f.den - self.den * f.num, self.den * f.den).reduce()

    def __mul__(self, f):
        return Fraction(self.num * f.num, self.den * f.den).reduce()

    def __truediv__(self, f):
        return Fraction(self.num * f.den, self.den * f.num).reduce()
    
def gramSchmidt(matrix):
    matrix = [[Fraction(x) for x in l] for l in matrix]
    gs = []
    for i, row in enumerate(matrix):
        gs.append(row)
        for j in range(i):
            gs[i] = [gs[i][k] - (dotp(gs[i], gs[j]) / dotp(gs[j], gs[j])) * gs[j][k] for k in range(len(gs[i]))]
    return gs
