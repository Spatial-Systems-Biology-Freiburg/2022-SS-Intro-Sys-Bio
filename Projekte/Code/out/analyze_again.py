import numpy as np
import json
from pathlib import Path
import multiprocessing as mp
import time
from collections import Counter
import matplotlib.pyplot as plt


def determine_if_pattern(y, component_index, xmax, ymax):    
    """
    Determines if a pattern was found
    
    Parameters
    ----------
    y : array-like
        result of the simulation at final time step with 3 dimensions
    component_indices : array-like
        contains all indices of components that are supposed to be analyzed for pattern formation
    xmax : int
        maximum domain of the simulation in x-direction
    ymax : int
        maximum domain of the simulation in y-direction
    
    Returns
    ----------
    s : bool
        If a pattern was not rejected, this parameter is True
    msg : string
        The message of how the parameter was (not) rejected
    
    ----------
    """
    # Determine if a pattern was created
    y_min = np.min(y, axis=(0,1))
    y_max = np.max(y, axis=(0,1))
    y_std = np.std(y, axis=(0,1))
    y_avg = np.average(y, axis=(0,1))
    
    # Check which nodes in the space count as peaks
    filt = (y > 0.5 * (y_max - y_min) + y_min)[:,:,component_index]
    n_peaks = np.sum(filt[:,:])

    # Reject if we only have fluctuations
    # if (y_std / (y_max - y_min))[component_index] > 1e-2:
    #     return False, "Rejected: y_std / (y_max - y_min) > 1e-2"
    # Reject if there are too many peaks
    if n_peaks >= y[:,:,component_index].size / 2.0:
        return False, "Rejected: n_peaks > xmax * ymax / 2"
    # Reject if number of peaks is not large enough
    if n_peaks < 4:
        return False, "Rejected: n_peaks < 4"
    # Reject if there are adjacent peaks in x or y direction
    if np.sum(filt[:-1]*filt[1:]) > 0:
        return False, "Rejected: Adjacent peaks"
    if np.sum(filt[:,:-1]*filt[:,1:]) > 0:
        return False, "Rejected: Adjacent peaks"
    # Otherwise assume that a pattern was successfully generated
    return True, "Not rejected"


def analyzer_function(filename):
    file = open(filename, "r")
    madata = json.load(file)

    solv_msg = madata["solving"]["message"] 
    lsa_succ = madata["LSA"]["succ"]        
    solv_succ = madata["solving"]["success"]

    # Now analyze final pattern again
    patt_succ, patt_msg = determine_if_pattern(
            np.array(madata["analysis"]["result_last"]),
            madata["analysis"]["component_index"],
            madata["model"]["xmax"],
            madata["model"]["xmax"]
    )
    file.close()
    return filename, lsa_succ, solv_succ, solv_msg, patt_succ, patt_msg


if __name__ == "__main__":
    start_time = time.time()
    files = sorted(Path().rglob("*.json"))
    print("Finding files:  {:06.2f}".format(time.time() - start_time))
    start_time = time.time()
    
    pool = mp.Pool()

    res = pool.starmap(
        analyzer_function,
        zip(
        files
        )
    )
    print("Analyzing files: {:06.2f}".format(time.time() - start_time))
    start_time = time.time()

    # Filter all results that have been reanalyzed and deemed true by new analysis
    new_files = [r[0] for r in res if r[4] == True]
 
    print("Analyzing files: {:06.2f}".format(time.time() - start_time))
    start_time = time.time()

    for filename in new_files:
        # Load the file and the last result
        file = open(filename, "r")
        data = json.load(file)
        
        print(filename)
        y = np.array(data["analysis"]["result_last"])[:,:,data["analysis"]["component_index"]]
        plt.imshow(y)
        plt.savefig("{}_im.png".format(filename))
        plt.clf()
