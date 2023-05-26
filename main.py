#!/usr/bin/env python3

import os
import argparse

from ase.io import read
import build_model
import runvasp
from build_model import Reaction


def execute_reaction(args):
    path, site, angle, height, reaction_num, runornot = args
    os.chdir(path)
    slab = read(os.path.join(path, 'slab/opt/CONTCAR'))
    reaction = build_model.Reaction(path=path,
                                    surf=slab,
                                    ads_points=slab[site].position,
                                    angle=angle,
                                    height=height)

    reaction_dict = {
        "oer": reaction.oer,
        "her": reaction.her,
        "hzor": reaction.hzor
    }

    if reaction_num in reaction_dict:
        reaction_dict[reaction_num]()
    else:
        raise ValueError("Wrong input reaction!")

    if runornot:
        runvasp.RunVasp(path=path, freq_index=build_model.freq_fix)()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reaction', type=str, default='oer')
    parser.add_argument('--path', type=str, default='.')
    parser.add_argument('--site', type=int, default=0)
    parser.add_argument('--angle', type=int, default=0)
    parser.add_argument('--height', type=float, default=1.7)
    parser.add_argument('--runornot', type=bool, default=True)

    args = parser.parse_args()
    execute_reaction([args.path, args.site, args.angle, args.height, args.reaction.lower(), args.runornot])
