# Catakit
The Python project is designed to provide an automated workflow for
simulating and analyzing various catalysis reactions,
including OER/ORR (Oxygen Evolution/Reduction Reaction),
HER (Hydrogen Evolution Reaction),
and HzOR (Horizon Zero Dawn of Reaction).  
The project enables users to build computational models,
run simulations, and generate step figures for visualizing the reaction mechanisms.

## Installing Dependencies

Before using this project, make sure you have the following libraries installed:

#### [pymatgen](https://pymatgen.org/installation.html)
#### [ASE](https://wiki.fysik.dtu.dk/ase/install.html)

## Usage
Git clone this project to some direction in your computer by using this command

```
cd DIR
git clone git@github.com:AlexHeustc/Catakit.git
```
`DIR`is your direction you want to save the project

***First you need to finish the slab structure optmization and put the files in direction like xxx/slab/opt/POSCAR POTCAR CONTCAR...***  

The VASP software is typically executed on a supercomputer, and a commonly used sbatch Bash script for running this project may look like this:
```
#!/bin/sh
#An example for MPI job.
#SBATCH -J job_name
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -p CPU
#SBATCH -N 2 -n 80
module load vasp/5.4.4
MPIRUN=mpirun #Intel mpi and Open MPI
python DIR/main.py --reaction hzor --path xxx/xxx --site 4 --angle 0 --runornot True --fixslab False
```
**Replace The `DIR/main.py` in this bash script as your real way to the project like `~/python_code/Catakit/main.py`**

**You can modify the [INCAR](./runvasp.py#44) and [the command to run vasp](runvasp.py#L59) to set your INCAR paramerter
to fit your computer, remember to modify the whole runvasp.py file!** 

**--reaction**  
This parameter is used to control the reaction you want to run. It accepts a string value and has a default_value of `hzor`. 
The optional arguments are `oer` `hzor` `her`.  

**--path**  
The path paremeter is the direction that you want to run the catalysis reaction  
The path is like `path/slab/opt`.  
The file tree willl be like :  
```angular2html
├── path/
│   ├── N2/
│   │   ├── freq/
│   │   ├── opt/
│   ├── N2H/
│   │   ├── freq/
│   │   ├── opt/
│   ├── N2H2/
│   │   ├── freq/
│   │   ├── opt/
│   ├── N2H3/
│   │   ├── freq/
│   │   ├── opt/
│   ├── N2H4/
│   │   ├── freq/
│   │   ├── opt/
│   ├── slab/
│   │   ├── opt/


```  
**--site(int)**  
Due to the python list index number starts from 0,the site argument is the real number-1,
like if you want to put the 5th atom to be the active center,it will be like `--site 4`  

**--angle(int,defult:0)**  
The angle you can rotate the adsorbates  

**--runornot(bool)**  
Whether to run Vasp or not

**--fixslab(bool)**  
Whether to fix slab when run vasp  

After run main.py you can copy the [step_figure.py](step_figure.py) to path,
and run the command `python step_figure.py`,then follow the guidance.  
The [step_figure.py](step_figure.py)  derived from modifying [jonyafei Python脚本提取数据绘制电化学台阶图](https://jonyafei.github.io/2021/06/20/python%E8%84%9A%E6%9C%AC%E6%8F%90%E5%8F%96%E6%95%B0%E6%8D%AE%E7%BB%98%E5%88%B6%E7%94%B5%E5%8C%96%E5%AD%A6%E5%8F%B0%E9%98%B6%E5%9B%BE/)
## Examples and Demos  
Here is an example:  
My bash file is:  
```angular2html
#!/bin/sh
#An example for MPI job.
#SBATCH -J job_name
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -p CPU-xxx
#SBATCH -N 3 -n 192

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.
MPIRUN=mpirun #Intel mpi and Open MPI
echo "Starting Time is `date`" >display
echo "Directory is `pwd`" >display
python ~/python_code/electrochemistry/main.py --reaction hzor --path /home/scms/xxx/xxx/FeNiRu --site 45 --angle 0 --runornot True --fixslab False
echo "Ending Time is `date`" >display
```
The file tree is like:
```angular2html
├── ./
│   ├── display
│   ├── job-19203.err
│   ├── job-19203.log
│   ├── job-19205.err
│   ├── job-19205.log
│   ├── mpi_job.sh
│   ├── N2/
│   │   ├── freq/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── DYNMAT
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── LOCPOT
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POT
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   ├── N2H/
│   │   ├── freq/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── DYNMAT
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── LOCPOT
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POT
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   ├── N2H2/
│   │   ├── freq/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── DYNMAT
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── LOCPOT
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POT
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   ├── N2H3/
│   │   ├── freq/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── DYNMAT
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── LOCPOT
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POT
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   ├── N2H4/
│   │   ├── freq/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── DYNMAT
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── LOCPOT
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POT
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
│   ├── slab/
│   │   ├── opt/
│   │   │   ├── CHG
│   │   │   ├── CHGCAR
│   │   │   ├── CONTCAR
│   │   │   ├── display
│   │   │   ├── DOSCAR
│   │   │   ├── EIGENVAL
│   │   │   ├── energy.txt
│   │   │   ├── IBZKPT
│   │   │   ├── INCAR
│   │   │   ├── job-17383.err
│   │   │   ├── job-17383.log
│   │   │   ├── KPOINTS
│   │   │   ├── mpi_job.sh
│   │   │   ├── OSZICAR
│   │   │   ├── OUTCAR
│   │   │   ├── PCDAT
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   ├── REPORT
│   │   │   ├── vasprun.xml
│   │   │   ├── WAVECAR
│   │   │   ├── XDATCAR
```
Then run step_figure.py and get
![这是图片](EnergyProfile.png "EnergyProfile")
