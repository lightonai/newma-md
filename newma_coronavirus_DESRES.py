# -*- coding: utf-8 -*-


"""
Script to recover the results for detecting the conformational changes of
SARS-CoV-2. Is it possible to choose
between real opu (to use it contact lighton at
 https://www.lighton.ai/contact-us/) and its synthetic version.
"""

import argparse

import numpy as np
import time

import onlinecp.algos as algos
import onlinecp.utils.feature_functions as feat
import onlinecp.utils.fastfood as ff

from MDAnalysis.lib.formats.libdcd import DCDFile

from lightonml.encoding.base import MultiThresholdEncoder, BinaryThresholdEncoder
from lightonml.projections.sklearn import OPUMap

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nmax', type=int,
                        default=7)
    #Trajectory studied: closed = trajectory1, opened = trajectory 2
    parser.add_argument('-trajclosed', type=str2bool, nargs='?',
                        const=True,default=True,
                        help='True for exploring traj (a), false for (b)')
    parser.add_argument('-repeat', type=int, default=0,
                        help='Augment the data to avoid artifacts at the '
                             'start')
    parser.add_argument('-full', type=str2bool, nargs='?',
                        const=True,default=False,
                        help='True for exploring all-inclusive traj')
    parser.add_argument('-ff', type=str2bool, nargs='?',
                        const=True, default=True,
                        help='True for performing NEWMA on CPU')
    args = parser.parse_args()

    nrepeat = args.repeat
    tc = args.trajclosed
    full = args.full
    newmaff = args.ff

    # If we study a trajectory including all atoms, a pyramid sketch is
    # needed on the OPU
    pyramid = full

    loadingtimestart = time.time()
    if tc:
        print('Exploring closed trajectory')
    else:
        print('Exploring open trajectory')

    # Path of the file containing the trajectories
    pathtraj = "/data/mldata/DESRES-Trajectory_"

    # Trajectory (a): initially closed
    nametraj1 = "sarscov2-10897136-no-water-no-ion-glueCA"

    # Trajectory (b): initially opened
    nametraj2 = "sarscov2-10897850-no-water-no-ion-glueCA"
    # Trajectory (b): all atoms
    namefull = "sarscov2-10897850-all-glueCA"
    pathtraj1 = pathtraj + nametraj1 + "/" + nametraj1 + "/" + nametraj1
    pathtraj2 = pathtraj + nametraj2 + "/" + nametraj2 + "/" + nametraj2
    pathfull = pathtraj + namefull + "/" + namefull + "/" + namefull

    # Initialize
    n0 = 0
    nmax = args.nmax

    if tc:
        path = pathtraj1
        print('Exploring closed trajectory')
    else:
        if full:
            path = pathfull
            print('Exploring opened trajectory - all')
        else:
            path = pathtraj2
            print('Exploring opened trajectory - no water no ion')

    traj0 = DCDFile(
        path + '-{:04d}'.format(n0) + ".dcd").readframes().xyz

    n_frames = traj0.shape[0]
    n_atoms = traj0.shape[1]
    n_feat = n_atoms * 3

    print(f"Number of time steps: {n_frames}. \nNumber of atoms: {n_atoms}." \
          f" \nFeatures: {n_feat}.")
    traj0 = traj0.reshape(n_frames, -1)

    traj = traj0

    if nrepeat > 0:
        print(f'Augmenting data {nrepeat} times.')
        for n in range(nrepeat):
            traj = np.concatenate((traj, traj0))

    for n in range(n0+1, nmax+1):
        print(n)
        pathfile = path + '-{:04d}'.format(n) + '.dcd'
        trajectory = DCDFile(pathfile).readframes().xyz
        trajectory = trajectory.reshape(n_frames, -1)
        traj = np.concatenate((traj, trajectory))


    tf = np.arange(1200 * ( 1 + n_frames * n0 ),
                   (nmax + 1) * 1200 * n_frames + 1, step=1200)

    print(f"Loading data took f{time.time() - loadingtimestart:.2f}s. ")

    # ----------
    # NEWMA RP OPU
    startingtime = time.time()

    # Pararameters for the NEWMA algorithm
    B = 250  # window size
    d = n_feat

    # Forget factors chosen with heuristic in the paper
    big_Lambda, small_lambda = algos.select_optimal_parameters(B)
    thres_ff = small_lambda
    # Number of random features is set automatically with this criterion
    m = int((1 / 4) / (small_lambda + big_Lambda) ** 2)
    m_OPU = 10 * m
    print(f'{m_OPU} random features.')


    if pyramid:
        # Pyramid sketch
        nc = 60000
        # Opening the OPU
        opu_mapping = OPUMap(n_components=nc)

        # First encoding
        encoder = BinaryThresholdEncoder(threshold_enc=10)
        minX, maxX = np.min(traj), np.max(traj)
        length_simulation = traj.shape[0]
        trajprojection = np.zeros((length_simulation, 3 * nc))

        for i in range(3):
            print(f"First sketch of the data... Step {i + 1} out of 3.")
            cut = n_atoms
            X = traj[:, cut * i:cut * (i + 1)]
            X = 255 * ((X - minX) / (maxX - minX))
            X = X.astype('uint8')
            Xencode = encoder.transform(X)
            del X

            X = opu_mapping.transform(Xencode)
            del Xencode
            trajprojection[:, i * nc: (i + 1) * nc] = X[:, :]
            del X

        print("Starting second sketch.")

        opu_mapping = OPUMap(n_components=m_OPU)
        n_levels = 38
        minX, maxX = np.min(trajprojection), np.max(trajprojection)
        print('Rescaling X')
        X = 255 * ((trajprojection - minX) / (maxX - minX))
        X = X.astype('uint8')

        print('Encoding X with the multi threshold encoder')
        thresholds = np.arange(n_levels, step=8)
        thresholds += 10
        encoderm = MultiThresholdEncoder(thresholds=thresholds)
        Xencode = encoderm.transform(X)
        del X

        X = opu_mapping.transform(Xencode)
        del Xencode
    else:
        # Regular sketch
        opu_mapping = OPUMap(n_components=m_OPU)
        n_levels = 38
        minX, maxX = np.min(traj), np.max(traj)
        print('Rescaling X')
        X = 255 * ((traj - minX) / (maxX - minX))
        X = X.astype('uint8')

        print('Encoding X with the multi threshold encoder')
        thresholds = np.arange(n_levels, step=8)
        thresholds += 10
        encoder = MultiThresholdEncoder(thresholds=thresholds)
        Xencode = encoder.transform(X)
        del X

        X = opu_mapping.transform(Xencode)
        del Xencode

    print(f"NEWMA OPU: Sketch took {time.time() - startingtime:.2f}s.")

    mult = 1.5
    detector = algos.NEWMA(X[0], forget_factor=big_Lambda,
                           feat_func=lambda x: x.astype('float32'),
                           forget_factor2=small_lambda,
                           adapt_forget_factor=thres_ff * mult,
                           thresholding_quantile=0.95,
                           dist_func=lambda z1, z2: np.linalg.norm(z1 - z2))
    detector.apply_to_data(X)

    n = 0
    detection_stat = np.array([i[0] for i in detector.stat_stored])[
                     int(10 * n):]  # padding
    online_th = np.array([i[1] for i in detector.stat_stored])[int(10 * n):]

    computation_duration = time.time() - startingtime
    print(f'NEWMA RP on OPU took {computation_duration:.2f}s.')

    # ----------
    # NEWMA FF on CPU
    if newmaff:
        startingtime = time.time()
        X = traj

        # Pararameters for the NEWMA algorithm
        B = 250
        d = n_feat

        # Forget factors chosen with heuristic in the paper
        big_Lambda, small_lambda = algos.select_optimal_parameters(B)
        thres_ff = small_lambda

        # number of random features is set automatically with this criterion
        m = 10 * int((1 / 4) / (small_lambda + big_Lambda) ** 2)

        choice_sigma = 'median'
        numel = 100
        data_sigma_estimate = X[:numel]

        W, sigmasq = feat.generate_frequencies(m, d, data=data_sigma_estimate,
                                               choice_sigma=choice_sigma)

        FF = ff.Fastfood(sigma=np.sqrt(sigmasq), n_components=m)
        FF.fit(X)
        X = FF.transform(X)

        detectorff = algos.NEWMA(X[0], forget_factor=big_Lambda,
                                 forget_factor2=small_lambda,
                                 adapt_forget_factor=0.5 * thres_ff)
        detectorff.apply_to_data(X)

        n = 0
        detection_statff = np.array([i[0] for i in detectorff.stat_stored])[
                           int(10 * n):]  # padding
        online_thff = np.array([i[1] for i in detectorff.stat_stored])[int(10 *
                                                                           n):]

        computation_duration = time.time() - startingtime
        print(f'NEWMA FF on CPU took {computation_duration:.2f}s.')


    detection_stat = detection_stat[nrepeat * n_frames:]
    online_th = online_th[nrepeat * n_frames:]

    tot_frame = n_frames * (nmax + nrepeat + 1)


    namefile = f"computetime_corona_totframe_{tot_frame}_natoms_{n_atoms}"

    # Save data for analysis in jupyter notebook.
    filename = f'newmaopu_totframe_{tot_frame}_natoms_{n_atoms}_closed_{tc}'
    np.savez_compressed(filename + f'_repeat{nrepeat}' + '.npz',
                        time=tf, detection_stat=detection_stat,
                        online_th=online_th)

