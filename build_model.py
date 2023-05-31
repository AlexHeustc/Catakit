import os
import shutil

import numpy as np
from ase import Atoms
from ase.build import sort
from ase.constraints import FixAtoms
from ase.io import read, write
from scipy.spatial import Delaunay

import runvasp

freq_fix = {}


class GetSite:
    def __init__(self, surf):
        self.surf = surf
        pos = self.surf.get_positions()
        tri = Delaunay(pos[:, :2])
        self.pos_nodes = pos[tri.simplices]

    def top(self):
        return self.surf.get_positions()

    def hollow(self):
        return np.mean(self.pos_nodes, axis=1)

    def bridge(self):
        bridge = []
        for triangle in self.pos_nodes:
            bridge.extend((triangle[i] + triangle[(i + 1) % 3]) / 2 for i in range(3))
        return np.array(bridge)


class Reaction:
    def __init__(self, path, surf, ads_points, angle, height=1.7, fixslab=1):
        self.surf = surf
        self.ads_points = ads_points
        self.height = np.array([0, 0, height])
        self.angle = angle
        self.fixslab = fixslab
        self.dir = []

    def build_model(self, adsorbates: dict):
        for label, positions in adsorbates.items():
            adsorbent = Atoms(label, positions + self.ads_points + self.height)
            # if isinstance(self.angle,int):
            adsorbent.rotate(self.angle, 'z', center=self.ads_points)
            # elif isinstance(self.angle,list) or isinstance(self.angle,tuple):
            # adsorbent.rotate(self.angle, center=self.ads_points)
            # adsorbent.rotate(self.angle[1], 'y', center=self.ads_points)
            # adsorbent.rotate(self.angle[2], 'z', center=self.ads_points)
            # 是否固定slab
            if self.fixslab != 0:
                self.surf.set_constraint(FixAtoms(indices=[atom.index for atom in self.surf]))
            atoms = self.surf + adsorbent

            atoms = sort(atoms)
            # 计算振动频率时固定
            freq_fix_indices = [atom.index for atom in atoms if atom.position.tolist() in self.surf.positions.tolist()]
            freq_fix[label] = freq_fix_indices

            directory = os.path.join(label, 'opt')
            os.makedirs(directory, exist_ok=True)
            write(os.path.join(directory, "POSCAR"), atoms)

    def her(self):
        ads = {'H': np.array([[0, 0, 0]])}
        height = 1.5
        self.build_model(ads)
        self.dir = ads.keys()

    def oer(self):
        ads = {'O': np.array([[0, 0, 0]]),
               'OH': np.array([[0, 0, 0], [0, 0, 0.9]]),
               'OOH': np.array([[0, 0, 0], [-0.1452, -0.8725, 1.18658], [-0.20479, -0.14158, 2.01980]])
               }
        self.build_model(ads)
        self.dir = ads.keys()

    def hzor(self):
        ads = {'N2H4': np.array([[0, 0, 0],
                                 [1.294395026, -0.005279826, 0.7034826],
                                 [-0.445146742, 0.884760161, 0.257637239],
                                 [-0.56252879, -0.725273111, 0.455652522],
                                 [1.831841066, 0.734513861, 0.234505592],
                                 [1.752386215, -0.867511408, 0.385781655]]
                                ),
               'N2H3': np.array([[0, 0, 0],
                                 [1.294395026, -0.005279826, 0.7034826],
                                 [-0.56252879, -0.725273111, 0.455652522],
                                 [1.831841066, 0.734513861, 0.234505592],
                                 [1.752386215, -0.867511408, 0.385781655]]),
               'N2H2': np.array([[0, 0, 0],
                                 [1.294395026, -0.005279826, 0.7034826],
                                 [-0.56252879, -0.725273111, 0.455652522],
                                 [1.831841066, 0.734513861, 0.234505592], ]),
               'N2H': np.array([[0, 0, 0],
                                [1.294395026, -0.005279826, 0.7034826],
                                [1.831841066, 0.734513861, 0.234505592], ]),
               'N2': np.array([[0, 0, 0],
                               [1.294395026, -0.005279826, 0.7034826], ])}

        self.build_model(ads)
        self.dir = ads.keys()

    def h2o(self):
        ads = {'OH2': np.array([[0, 0, 0], [0, -0.79068950, 0.61221720], [0, 0.79068950, 0.61221720]])}
        self.build_model(ads)
        self.dir = ads.keys()

    def o2(self):
        ads = {'O2': np.array([[0, 0, 0], [0, 0, 1.23]])}
        self.build_model(ads)
        self.dir = ads.keys()

    def iocn(self, iocn_symbol: str):
        ads = {iocn_symbol: np.array([[0, 0, 0]])}
        height = 1.5
        self.build_model(ads)
        self.dir = ads.keys()
