import numpy as np
import time
import matplotlib.pyplot as plt
import itertools as it
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
import multiprocessing as mp
import os


def __save_result_plot(i, u, index, min, max, start_time=None, output_folder=Path("out")):
    if start_time!=None:
        print("[{: >8.4f}s] Saving Plots ...".format(time.time()-start_time), end="\r")
    fig, ax = plt.subplots()
    im = ax.imshow(
        u,
        vmin=min[index],
        vmax=max[index],
        cmap='YlGn',
        interpolation='spline36'
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    fig.savefig(output_folder / "image_{:08d}.png".format(i))
    plt.close(fig)


def save_plots(res, component_index, step, threads=os.cpu_count(), output_folder=Path("./out/")):
    start_time = time.time()
    with mp.Pool(threads) as p:
        p.starmap(__save_result_plot, zip(
            range(0, res.shape[3], step),
            [res[:,:,component_index,i] for i in range(0, res.shape[3], step)],
            it.repeat(component_index),
            it.repeat(np.min(res, axis=(0, 2, 3))),
            it.repeat(np.max(res, axis=(0, 2, 3))),
            it.repeat(start_time),
            it.repeat(output_folder)
        ))
    print("[{: >8.4f}s] Saved all Plots".format(time.time()-start_time))