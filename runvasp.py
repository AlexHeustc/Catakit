import os
import shutil
import subprocess
import warnings

from ase.constraints import FixAtoms
from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.ase import AseAtomsAdaptor

warnings.filterwarnings('ignore')


def check_vasp_status(dir: str) -> bool:
    # 检查OUTCAR文件的最后几行是否包含已经收敛的关键词
    keyword_list = ["reached required accuracy"]
    # first judge whether the file exists
    if os.path.exists(os.path.join(dir, 'OUTCAR')):
        with open(f"{dir}/OUTCAR", "r") as f:
            lines = f.readlines()[-10:]  # 读取OUTCAR文件的最后10行
            for line in lines:
                for keyword in keyword_list:
                    if keyword in line:
                        return True
                    else:
                        return False
    else:
        return False


class RunVasp:
    def __init__(self, path, freq_index: dict):
        self.path = path
        self.dirs = [os.path.join(path, dir_key) for dir_key in freq_index.keys()]
        self.freq_index = freq_index

    def __call__(self, *args, **kwargs):
        self.perform_optimization()
        self.calculate_frequency()

    def perform_optimization(self):
        custom_settings = {
            "NELMIN": 5, "NELM": 300, "NCORE": 8,
            "EDIFF": 1e-4, "EDIFFG": -0.05,
            "ISMEAR": 0, "SIGMA": 0.05,
            "NSW": 500, "PREC": "Normal", "ISIF": 2, "ISPIN": 2, "LORBIT": None, "AlGO": "Fast", "LCHARG": False,
            "LWAVE": False, "IVDW": 11,
        }
        for dir in self.dirs:
            # check the optimization
            opt_path = os.path.join(dir, 'opt')
            if check_vasp_status(opt_path):
                continue
            else:
                os.chdir(opt_path)
                relax = MPRelaxSet(structure=Structure.from_file('POSCAR'), user_incar_settings=custom_settings)
                relax.write_input('.')
                os.system("mpirun vasp_std >>display")
                os.chdir(self.path)

    def calculate_frequency(self):
        custom_settings = {
            "NELMIN": 5, "NELM": 300, "NCORE": 8, "ISMEAR": 0, "SIGMA": 0.05, "NSW": 1, "PREC": "Normal",
            "EDIFF": 1e-4, "EDIFFG": -0.05,
            "ISIF": 2, "ISPIN": 2, "ISYM": 0, "LREAL": "A", "LORBIT": None,
            "AlGO": "Fast", "LCHARG": False, "LWAVE": "F", "IVDW": 11, "NFREE": 2, "POTIM": 0.015,
            "IBRION": 5, "LAECHG": False}

        for dir in self.dirs:
            freq_path = os.path.join(dir, 'freq')
            # 进入工作目录
            os.chdir(dir)
            freq_path = os.path.join(dir, 'freq')
            os.makedirs(freq_path, exist_ok=True)
            if check_vasp_status(freq_path):
                continue
            else:
                os.chdir(freq_path)
                os.system("cp ../opt/CONTCAR POSCAR")
                atoms = read("POSCAR")
                c = FixAtoms(indices=self.freq_index[f"{dir.split('/')[-1]}"])
                atoms.set_constraint(constraint=c)

                zep = MPStaticSet(structure=AseAtomsAdaptor.get_structure(atoms),
                                  user_incar_settings=custom_settings)
                zep.write_input('.')
                os.system("mpirun vasp_std >>display")

                os.chdir(self.path)
