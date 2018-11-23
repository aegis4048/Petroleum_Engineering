class Variogram2D():

    def __init__(self, foo, bar, azi):
        self.a = []
        self.b = []

        for i in range(2):
            c = np.array(foo)
            d = np.array(bar)

            c.angle = azi[i]
            d.angle = azi[i]

            self.a.append(c)
            self.b.append(d)

