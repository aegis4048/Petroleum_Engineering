class Values2D:
    def __init__(self, values, azimuth):
        self.values = values
        self.azimuth = azimuth


class Variogram2D:

    def __init__(self, x, y, z, dist_max, azi, azi_tol, lag_dist, lag_tol):

        assert len(x) == len(y) == len(z), "x, y, and z must have the same dimension"
        assert 0 <= azi <= 360, "azimuth must have a range of [0, 360]"
        assert 0 <= azi_tol <= 90, "azimuth tolerance must have a range of [0, 90]"

        if azi == 90:
            self.isotropy = True
        else:
            self.isotropy = False

        self.azi_tol = azi_tol
        self.lag_dist = lag_dist
        self.lag_tol = lag_tol

        self.npairs_2D = []
        self.sumsq_2D = []
        self.gamma_2D = []

        # number of data points
        self.n = len(x)

        # lag distances on which variogram gamma values will be computed
        self.lags = np.array([i * 5 for i in range(int(dist_max / lag_dist) + 1)])

        # 1-D distances between two points on a coordinate system
        self.dx = np.array([[x[i] - x[j] for j in range(self.n)] for i in range(self.n)])
        self.dy = np.array([[y[i] - y[j] for j in range(self.n)] for i in range(self.n)])

        # 2-D distances between two data points on a coordinate system
        self.distances = np.array(
            [[np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) for j in range(self.n)] for i in range(self.n)])

        # difference in angle between two data points (degrees)
        self.dazimuth = np.array([[self.calc_theta(i, j) for j in range(self.n)] for i in range(self.n)])

        # squared differences between two data points values at a certain lag distance
        self.pairwise_sq_diff = np.array([[(z[i] - z[j]) ** 2 for j in range(self.n)] for i in range(self.n)])

        # user defined azimuth in two directions. In case of isotropy, azimuth is the same for both directions
        if self.isotropy:
            self.azi = [self.ensure_positive_azimuth(azi), self.ensure_positive_azimuth(azi)]
        else:
            self.azi = [self.ensure_positive_azimuth(azi), self.ensure_positive_azimuth(azi + 90)]

        # variogram gamma values in two directions. They are the same for both directions in case of isotropy
        npairs = Values2D
        sumsq = Values2D
        gamma = Values2D
        for i in range(2):
            npairs.values = np.array([self.calc_npairs(lag, self.azi[i]) for lag in self.lags])
            sumsq.values = np.array([self.calc_sumsq(lag, self.azi[i]) for lag in self.lags])
            gamma.values = np.divide(sumsq.values, npairs.values * 2)

            npairs.azimuth = self.azi[i]
            sumsq.azimuth = self.azi[i]
            gamma.azimuth = self.azi[i]

            self.npairs_2D.append(npairs)
            self.sumsq_2D.append(sumsq)
            self.gamma_2D.append(gamma)

    # ensures positive azimuth between two data points.
    # In variograms, positive y-axis of a 2-D plane is 0 degree azimuth. Azimuth increases clockwise
    def ensure_positive_azimuth(self, azi):
        if azi + self.azi_tol > 360:
            return azi - 180
        elif azi - self.azi_tol < 0:
            return azi + 180
        else:
            return azi

    def calc_theta(self, i, j):
        if self.dx[i, j] > 0:
            theta = np.degrees(np.pi / 2 - np.arctan(self.dy[i, j] / self.dx[i, j]))
        elif self.dx[i, j] < 0:
            theta = np.degrees(np.pi * 1.5 - np.arctan(self.dy[i, j] / self.dx[i, j]))
        else:
            if self.dy[i, j] > 0:
                theta = 0
            elif self.dy[i, j] < 0:
                theta = 180
            else:
                theta = 0
        return theta

    def is_within_azi_lag_tolerance(self, data_azi, data_lag, lag, azi):
        is_azi_tolerated = azi - self.azi_tol <= data_azi <= azi + self.azi_tol
        is_lag_tolerated = lag - self.lag_tol <= data_lag <= lag + self.lag_tol
        return is_azi_tolerated and is_lag_tolerated

    def calc_npairs(self, lag, azi):
        npairs = 0
        for i in range(self.n):
            for j in range(self.n):
                if self.is_within_azi_lag_tolerance(self.dazimuth[i, j], self.distances[i, j], lag, azi):
                    npairs += 1
        return npairs

    def calc_sumsq(self, lag, azi):
        sumsq = 0
        for i in range(self.n):
            for j in range(self.n):
                if self.is_within_azi_lag_tolerance(self.dazimuth[i, j], self.distances[i, j], lag, azi):
                    sumsq += self.pairwise_sq_diff[i, j]
        return sumsq


temp = Variogram2D(x, y, z, 100, 340, 45, 5, 10)