import pygmsh
import meshio
import numpy as np
from CoolProp.CoolProp import PropsSI
import pandas as pd

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
    def __init__(self, T_min, T_max, D_min, fluid, nT, nD):
        self.T_min = T_min
        self.T_max = T_max
        self.D_min = D_min
        self.D_max = PropsSI('D',
                             "T", PropsSI("Tcrit", fluid),
                             'P', PropsSI("Pcrit", fluid), fluid)
        self.fluid = fluid
        self.nT = nT
        self.nD = nD

    def create_table(self,req_prop):
        self._create_mesh()
        self._generate_point_data(req_prop)
        self._convert_mesh_to_linear()

    def write_table(self, filename):
        meshio.write(filename, self.mesh, write_binary=False)

    def _create_mesh(self):
        self.geom = pygmsh.built_in.Geometry()
        c1, c2, c3, c4 = self._get_corners()
        b, r, t, l = self._get_lines(c1, c2, c3, c4)
        s = self._get_surface([b,r,t,l])
        self.mesh = self._get_mesh(s)

    def _generate_point_data(self,req_prop):
        self.mesh.point_data = self._get_point_data(req_prop)

    def _get_corners(self):
        c1 = self.geom.add_point(np.array([np.log10(self.D_min), self.T_min, 0.0]))
        c2 = self.geom.add_point(np.array([np.log10(self.D_max), self.T_min, 0.0]))
        c3 = self.geom.add_point(np.array([np.log10(self.D_max), self.T_max, 0.0]))
        c4 = self.geom.add_point(np.array([np.log10(self.D_min), self.T_max, 0.0]))
        return c1, c2, c3, c4


    def _get_lines(self, c1, c2, c3, c4):
        b = self.geom.add_line(c1, c2)
        r = self.geom.add_line(c2, c3)
        t = self.geom.add_line(c3, c4)
        l = self.geom.add_line(c4, c1)
        return b, r, t, l

    def _get_surface(self, lines):
        ll = self.geom.add_line_loop(lines)
        s = self.geom.add_surface(ll)
        return s

    def _get_mesh(self,s):
        self.geom.set_transfinite_surface(s, size=[self.nD, self.nT])
        self.geom.add_raw_code('Recombine Surface {%s};' % s.id)
        self.geom.add_raw_code("Mesh.Algorithm = 9;")
        return pygmsh.generate_mesh(self.geom, dim=2)

    def _convert_mesh_to_linear(self):
        for pt in self.mesh.points:
            pt[0] = np.power(10, pt[0])

    def _convert_df_to_dict(self,df):
        d = dict()
        for col in df.columns:
            d[col]=df[col].values
        return d

    def _get_point_data(self,req_prop):
        df = pd.DataFrame(self.mesh.points, columns=['logD', 'T', 'dummy'])
        df['tpoints']=df.apply(lambda x: ThermodynamicPoint(np.power(10,x['logD']),x['T'],self.fluid),axis=1)
        for prop in req_prop:
            df[prop]=df['tpoints'].apply(lambda x: getattr(x,prop))
        return self._convert_df_to_dict(df[req_prop])

if __name__=="__main__":
    filename = "table.vtk"
    req_prop = ['H', 'T', 'P', 'D', 'U', 'V','A', 'Q']
    fluid = 'Toluene'
    T_min = 400
    T_max = 650
    D_min = 0.001

    table = Table(T_min=T_min, T_max=T_max, fluid=fluid,nD=100,nT=10, D_min=D_min)
    table.create_table(req_prop)
    table.write_table(filename)
