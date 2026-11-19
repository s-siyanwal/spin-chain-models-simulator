class Moments4:
    def __init__(self): self.n=0; self.m1=self.m2=self.m3=self.m4=0.0
    def update(self,x):
        n1=self.n; self.n+=1; d=x-self.m1; dn=d/self.n; dn2=dn*dn; term=d*dn*n1
        self.m4 += term*dn2*(self.n*self.n-3*self.n+3)+6*dn2*self.m2-4*dn*self.m3
        self.m3 += term*dn*(self.n-2)-3*dn*self.m2
        self.m2 += term; self.m1 += dn
    def add_symmetrized(self,s): self.update(abs(s)); self.update(-abs(s))
    @property
    def kurtosis(self): return self.n*self.m4/(self.m2*self.m2) if self.n>3 and self.m2 else float('nan')
