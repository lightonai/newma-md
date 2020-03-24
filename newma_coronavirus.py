# -*- coding: utf-8 -*-


"""
Script to recover the results for detecting the conformational changes of
SARS-CoV-2. Is it possible to choose
between real opu (to use it contact lighton at
 https://www.lighton.ai/contact-us/) and its synthetic version.
"""


import numpy as np
import time

import onlinecp.algos as algos
import onlinecp.utils.feature_functions as feat
import onlinecp.utils.fastfood as ff

from MDAnalysis import Universe

from lightonml.encoding.base import MultiThresholdEncoder
from lightonml.projections.sklearn import OPUMap


if __name__ == '__main__':
    # Path of the file containing the trajectory
    pathtraj = "./11764158/QHD_LP1/QHD_LP1/"

    # Importing the trajectory using MDAnalysis
    # Note: different tools might be needed depending on the format of the data
    universe = Universe(pathtraj + 'QHD_LP1.gro', pathtraj + 'QHD_LP1.xtc')
    all_atoms = universe.select_atoms("all")

    traj = np.asarray([all_atoms.positions for frame in universe.trajectory],
                        dtype='f8')
    traj = traj.reshape(traj.shape[0], -1)

    timestart = 0
    timestop = traj.shape[0]

    n_feat = traj.shape[1]
    n_atoms = np.int(n_feat/3)

    print(f"Number of atoms: {n_atoms}.")
    print(f"Number of features: {n_features}.")
    print(f"Number of timeframes: {timestop - timestart}")

    # ----------
    # NEWMA RP OPU
    startingtime = time.time()

    # Pararameters for the NEWMA algorithm
    B = 80  # window size
    d = n_feat

    # Forget factors chosen with heuristic in the paper
    big_Lambda, small_lambda = algos.select_optimal_parameters(B)
    thres_ff = small_lambda
    # Number of random features is set automatically with this criterion
    m = int((1 / 4) / (small_lambda + big_Lambda) ** 2)
    m_OPU = 10 * m
    print(f'{m_OPU} random features.')

    # Opening the OPU
    opu_mapping = OPUMap(n_components=m_OPU)

    # Encoding the data
    n_levels = 38
    minX, maxX = np.min(traj), np.max(traj)
    print('Rescaling X')
    X = 255 * ((traj - minX) / (maxX - minX))
    X = X.astype('uint8')

    print('Encoding X with the multi threshold encoder')
    thresholds = np.arange(n_levels, step=5)
    thresholds += 10
    encoder = MultiThresholdEncoder(thresholds=thresholds)
    Xencode = encoder.transform(X)
    del X

    X = opu_mapping.transform(Xencode)
    del Xencode

    # NEWMA
    mult = 0.5
    detector = algos.NEWMA(X[0], forget_factor=big_Lambda,
                            feat_func=lambda x: x.astype('float32'),
                            forget_factor2=small_lambda,
                            adapt_forget_factor=thres_ff*mult,
                            thresholding_quantile=0.95,
                            dist_func=lambda z1, z2: np.linalg.norm(z1 - z2))
    detector.apply_to_data(X)

    detection_stat = np.array([i[0] for i in detector.stat_stored])
    online_th = np.array([i[1] for i in detector.stat_stored])

    computation_duration = time.time() - startingtime
    print(f'Computation RP OPU took {computation_duration:.2f}s.')

    # ----------
    # NEWMA RFF CPU
    startingtime = time.time()
    X = traj

    # Pararameters for the NEWMA algorithm
    B = 80  # window size
    d = n_feat

    # Forget factors chosen with heuristic in the paper
    big_Lambda, small_lambda = algos.select_optimal_parameters(B)
    thres_ff = small_lambda
    # number of random features is set automatically with this criterion
    m = 10 * int((1 / 4) / (small_lambda + big_Lambda) ** 2)

    choice_sigma = 'median'
    numel = 100
    data_sigma_estimate = X[:numel]  # data for median trick to estimate sigma

    W, sigmasq = feat.generate_frequencies(m, d, data=data_sigma_estimate,
                                           choice_sigma=choice_sigma)

    def feat_func(x):
        return feat.fourier_feat(x, W)

    # NEWMA
    detectorb = algos.NEWMA(X[0], forget_factor=big_Lambda,
                            forget_factor2=small_lambda,
                            feat_func=feat_func,
                            adapt_forget_factor=0.5 * thres_ff)
    detectorb.apply_to_data(X)

    detection_statrff = np.array([i[0] for i in detectorb.stat_stored])
    online_thrff = np.array([i[1] for i in detectorb.stat_stored])

    computation_duration = time.time() - startingtime
    print(f'Computation RFF took {computation_duration:.2f}s.')

    # ---------
    # NEWMA FF CPU
    startingtime = time.time()
    X = traj

    # Parameters for the NEWMA algorithm
    B = 80 # window size
    d = n_feat

    # Forget factors chosen with heuristic in the paper
    big_Lambda, small_lambda = algos.select_optimal_parameters(B)
    thres_ff = small_lambda
    # number of random features is set automatically with this criterion
    m = 10 * int((1 / 4) / (small_lambda + big_Lambda) ** 2)

    choice_sigma = 'median'
    numel = 100
    data_sigma_estimate = X[:numel] # data for median trick to estimate sigma

    W, sigmasq = feat.generate_frequencies(m, d, data=data_sigma_estimate,
                                           choice_sigma=choice_sigma)

    FF = ff.Fastfood(sigma=np.sqrt(sigmasq), n_components=m)
    FF.fit(X)
    X = FF.transform(X)

    # NEWMA
    detectorff = algos.NEWMA(X[0], forget_factor=big_Lambda,
                            forget_factor2=small_lambda,
                            adapt_forget_factor=0.5 * thres_ff)
    detectorff.apply_to_data(X)

    detection_statff = np.array([i[0] for i in detectorff.stat_stored])
    online_thff = np.array([i[1] for i in detectorff.stat_stored])

    computation_duration = time.time() - startingtime
    print(f'Computation FF took {computation_duration:.2f}s.')



    # Save data for analysis in jupyter notebook.
    np.savez_compressed(f'newmaopu_corona_natoms_{n_atoms}.npz',
                        detection_stat=detection_stat,
                        online_th=online_th)

