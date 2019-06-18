import pygmsh
import meshio
import numpy as np
from CoolProp.CoolProp import PropsSI
import pandas as pd
import matplotlib.pyplot as plt

class ThermodynamicPoint(object):
    def __init__(self, D, T, fluid):
        self.D = D
        self.T = T
        self.P = PropsSI('P','D',self.D,'T',self.T,fluid)
        self.H = PropsSI('H', 'D', self.D, 'T', self.T,fluid)
        self.V = PropsSI('V', 'D', self.D, 'T', self.T,fluid)
        self.U = PropsSI('U', 'D', self.D, 'T', self.T,fluid)
        self.Q = PropsSI('Q', 'D', self.D, 'T', self.T,fluid)
        if self.Q < 0.0:
            self.Q = 1.0
        if self.Q == 1.0:
            self.A = PropsSI('A', 'D', self.D, 'T', self.T, fluid)
        else:
            self.A = PropsSI('A', 'Q', 1.0, 'T', self.T,fluid)*(self.Q) + PropsSI('A', 'Q', 0.0, 'T', self.T,fluid)*(1-self.Q)

class Table(object):
    def __init__(self, T_min, T_max, fluid,
                 size_Crit,
                 size_Tmax_Dmax,
                 size_Tmax_Dmin,
                 size_Tmin_Dmin):
        self.T_min = T_min
        self.T_max = T_max
        self.D_min = PropsSI("D", "Q", 1.0,"T", T_min, fluid)
        self.D_max = PropsSI('D',
                             "T", PropsSI("Tcrit", fluid),
                             'P', PropsSI("Pcrit", fluid)*0.99, fluid)

        self.fluid = fluid
        self.size_Crit = size_Crit
        self.size_Tmax_Dmax = size_Tmax_Dmax
        self.size_Tmax_Dmin = size_Tmax_Dmin
        self.size_Tmin_Dmin = size_Tmin_Dmin


    def create_table(self,req_prop):
        self._create_mesh()
        self._generate_point_data(req_prop)
        # self._convert_mesh_to_linear()

    def write_table(self, filename):
        meshio.write(filename, self.mesh, write_binary=False)


    def _create_boundaries(self):
        n = 100
        D = np.geomspace(self.D_max, self.D_min, n)
        bot = pd.DataFrame(data=D, columns=['D'])
        bot['T'] = bot['D'].apply(lambda x: PropsSI("T", "D", x, "Q", 1.0, fluid))
        bot = self._add_nondim(bot)
        bot['size'] = bot.apply(lambda x: self.size_Crit +
                                          (self.size_Tmin_Dmin - self.size_Crit) * (1.0 - x['Dstar']), axis=1)

        # top
        top = pd.DataFrame(data=D, columns=['D'])
        top['T'] = self.T_max
        top = self._add_nondim(top)
        top['size'] = top.apply(lambda x: self.size_Tmax_Dmax +
                                          (self.size_Tmax_Dmin - self.size_Tmax_Dmax) * (1.0 - x['Dstar']), axis=1)

        # left
        left = pd.DataFrame(data=np.linspace(self.T_min, self.T_max, n), columns=['T'])
        left['D'] = self.D_min
        left = self._add_nondim(left)
        left['size'] = left.apply(lambda x: self.size_Tmax_Dmin +
                                            (self.size_Tmin_Dmin - self.size_Tmax_Dmin) * (1 - x['Tstar']), axis=1)

        ##right
        right = pd.DataFrame(data=np.linspace(self.T_max, bot['T'].iloc[0], n), columns=['T'])
        right['D'] = self.D_max
        right = self._add_nondim(right)
        right['size'] = right.apply(lambda x: self.size_Crit + (self.size_Tmax_Dmax - self.size_Crit) *
                                                               (x['Tstar'] - bot['Tstar'].iloc[0]), axis=1)


        return bot, top, right, left

    def _create_mesh(self):
        self.geom = pygmsh.built_in.Geometry()
        self.bot, self.top, self.right, self.left = self._create_boundaries()
        self.right = self._add_points(self.right)
        self.left = self._add_points(self.left)
        self.top = self._add_points(self.top)
        self.bot = self._add_points(self.bot)
        b = self.geom.add_spline(self.bot['points'].values)
        l = self.geom.add_line(self.bot['points'].iloc[-1], self.   left['points'].iloc[-1])
        t = self.geom.add_line(self.left['points'].iloc[-1], self.right['points'].iloc[0])
        r = self.geom.add_line(self.right['points'].iloc[0], self.bot['points'].iloc[0])
        ll = self.geom.add_line_loop([r, b, l, t])
        s = self.geom.add_surface(ll)
        self.mesh = pygmsh.generate_mesh(self.geom, dim=2)

    def _generate_point_data(self,req_prop):
        self.mesh.point_data = self._get_point_data(req_prop)


    def _get_point_data(self,req_prop):
        df = pd.DataFrame(self.mesh.points, columns=['Dstar', 'Tstar', 'dummy'])
        df['logD'] = df['Dstar'].apply(lambda x: x * (np.log10(self.D_max) - np.log10(self.D_min)) + np.log10(self.D_min))
        df['T'] = df['Tstar'].apply(lambda x: (x * (self.T_max - self.T_min) + self.T_min))
        df['tpoints'] = df.apply(lambda x: ThermodynamicPoint(np.power(10, x['logD']), x['T'], fluid), axis=1)
        for prop in req_prop:
            df[prop]=df['tpoints'].apply(lambda x: getattr(x,prop))
        return self._convert_df_to_dict(df[req_prop])
    
    def _convert_mesh_to_linear(self):
        for pt in self.mesh.points:
            pt[0] = np.power(10, pt[0])

    def _convert_df_to_dict(self,df):
        d = dict()
        for col in df.columns:
            d[col]=df[col].values
        return d

    def _add_nondim(self,df):
        df['Tstar']=df['T'].apply(lambda x: (x - self.T_min)/(self.T_max-self.T_min))
        df['Dstar'] = df['D'].apply(lambda x: (np.log10(x) - np.log10(self.D_min)) / (np.log10(self.D_max) - np.log10(self.D_min)))
        return df

    def _add_points(self,df):
        df['points']=df.apply(lambda x: self.geom.add_point(np.array([x['Dstar'],x['Tstar'], 0.0]),x['size']),axis=1)
        return df

if __name__=="__main__":
    filename = "table.vtk"


    #inputs
    fluid = 'Toluene'
    T_min = 300
    T_max = 650
    req_prop = ['H', 'T', 'P', 'D', 'U', 'V', 'A', 'Q']

    #mesh inputs
    size_Tmax_Dmax = 0.01
    size_Crit = 0.001
    size_Tmax_Dmin = 0.04
    size_Tmin_Dmin = 0.02

    #create unstructured table
    table = Table(T_min, T_max, fluid,
                 size_Crit,
                 size_Tmax_Dmax,
                 size_Tmax_Dmin,
                 size_Tmin_Dmin)

    table.create_table(req_prop)
    table.write_table("table.vtk")


    #plotting
    fig, ax = plt.subplots()
    table.bot.plot(x='D',y='T', logx=True, ax=ax,marker='x', label='bot_bound')
    table.left.plot(x='D', y='T', ax=ax, marker='x', label='left_bound')
    table.top.plot(x='D', y='T', ax=ax, marker='x', label='top_bound')
    table.right.plot(x='D', y='T', ax=ax, marker='x', label='right_bound')
    ax.set_xlim([table.D_min*0.95, table.D_max*1.05])
    ax.set_xlabel("Density*")
    ax.set_ylabel("Temperature*")

    fig, ax = plt.subplots()
    table.bot.plot(x='Dstar', y='Tstar', logx=False, ax=ax, marker='x', label='bot_bound')
    table.left.plot(x='Dstar', y='Tstar', ax=ax, marker='x', label='left_bound')
    table.top.plot(x='Dstar', y='Tstar', ax=ax, marker='x', label='top_bound')
    table.right.plot(x='Dstar', y='Tstar', ax=ax, marker='x', label='right_bound')
    ax.set_xlim([-0.1, 1.1])
    ax.set_xlabel("D^*")
    ax.set_ylabel("T^*")

    plt.show()