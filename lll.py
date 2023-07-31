from math import gcd, floor

def dotp(a, b): return sum([x * y for x, y in zip(a,b)], start=Fraction(0))

class Fraction():
    def __init__(self, num, den=1):
        self.num = num
        self.den = den
    
    def reduce(self):
        g = gcd(self.num, self.den)
        if g != 1:
            self.num //= g
            self.den //= g
        return self
    
    def near(self):
        return Fraction(floor(self.num/self.den + 0.5))


    def __repr__(self): return f"({self.num}/{self.den})"

    def __add__(self, f):
        return Fraction(self.num * f.den + self.den * f.num, self.den * f.den).reduce()

    def __sub__(self, f):
        return Fraction(self.num * f.den - self.den * f.num, self.den * f.den).reduce()

    def __mul__(self, f):
        return Fraction(self.num * f.num, self.den * f.den).reduce()

    def __truediv__(self, f):
        return Fraction(self.num * f.den, self.den * f.num).reduce()

    def __pow__(self, n):
        return Fraction(self.num**n, self.den**n).reduce()

    def __lt__(self, f):
        return (self.num * f.den < self.den * f.num)
     
    def __gt__(self, f):
        return (self.num * f.den > self.den * f.num)

    def __eq__(self, f):
        return (self.num * f.den == self.den * f.num)
    
def gramSchmidt(matrix):
    gs = []
    for i, row in enumerate(matrix):
        gs.append(row)
        for j in range(i):
            gs[i] = [gs[i][k] - (dotp(gs[i], gs[j]) / dotp(gs[j], gs[j])) * gs[j][k] for k in range(len(gs[i]))]
    return gs

def lll(matrix, delta):
    N = len(matrix)

    b = [[Fraction(x) for x in l] for l in matrix]
    b_ = gramSchmidt(b)
    print("before", b)
    print("gso", b_)

    def mu(i, j): return dotp(b[i], b_[j]) / dotp(b_[j], b_[j])

    k = 1
    while (k < N):
        for j in reversed(range(k)):
            print("Iter", k, j)
            _ = input()
            m = mu(k, j)
            print("mu", m)
            _ = input()
            if m > Fraction(1,2):
                b[k] = [b[k][i] - m.near() * b[j][i] for i in range(N)]
                print(f"changed row {k} to {b[k]}")
                _ = input()
                b_ = gramSchmidt(b)
                print("updated gso", b_)
                _ = input()
        if (dotp(b_[k], b_[k]) > (delta - mu(k, k-1)**2) * dotp(b_[k-1], b_[k-1])):
            print("incrementing k")
            _ = input()
            k += 1
        else:
            b[k], b[k-1] = b[k-1], b[k]
            print("swapped")
            _ = input()
            b_ = gramSchmidt(b)
            print("updated gso", b_)
            _ = input()
            k = max(k-1, 1)
    return b
