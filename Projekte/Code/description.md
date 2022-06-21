# PDE Solver
## Requirements
In order to use the provided pde solver, the user must have working installations of the following python packages.
```python
import numpy
import scipy
import matplotlib
```
This packages was developed with python `3.10.5`, numpy `1.22.4` and scipy `1.8.1`.
However, we expect most versions from python `3.5` upwards to work.
## Usage
We first need to import needed modules.
```python
from PDE_Solver import PDE_Solver, save_result_plot, save_plots
import numpy as np
```
Afterwards we define the reaction kinetics governing our system
```python
def ODE(t, y, D, ind, k):
    dydt = np.zeros(y.shape)
    A = y[ind+0]
    B = y[ind+1]
    dydt[ind]   = k[0] - k[1]*A + k[2]*A**2*B   + k[3]*np.dot(D,A)
    dydt[ind+1] = k[4] - k[2]*A**2*B            + k[5]*np.dot(D,B)
    return dydt
```
To use the pde solver, we need to define grid size, #components, lower and upper bound for time interval, number of time steps, initial values, type of boundary conditions, grid-type and kinetics parameters
```python
xmax=20
ymax=20
NVar=2
t_span = (0, 100)
t_num = 500
y0 = np.random.uniform(low=600, high=1000, size=xmax*ymax*NVar)
bndcondition="zeroflux"
celltype="quadratic"

k = [
	10.0, 0.1, 4.83e-07, 100.0/40**2,
	80.0, 5000.0/40**2
]
```
Now we need some helper variables such as time points to solve for, coupling matrix and index generator.
```python
t_eval = np.linspace(t_span[0], t_span[1], t_num)
D=couplingMatrix(xmax, ymax, bndcondition, celltype)
ind = IJKth(1, np.arange(ymax), np.arange(xmax), ymax, NVar)
```
Now we can solve the equations
```python
start_time = time.time()
print("[{: >8.4f}] Solving ...".format(0), end="\r")
sol = solve_ivp(
	lambda t, y: ODE(t, y, D, ind, k),
	t_span,
	y0,
	method='Radau',
	jac_sparsity=jpat(D,ind),
	vectorized=True,
	t_eval=t_eval
)
res = sol.y.reshape((xmax, ymax, NVar, len(t_eval)))
print("[{: >8.4f}s] Solving Done".format(time.time()-start_time))
```
We can save plots by using the predefined method(s)
```python
# This will save the singular result at time step i of the results np.array
# The argument start_time can be neglected when calling individually
save_result_plot(i, res, index, min, max, start_time=None, output_folder=Path("./out/"))

# This will save every n_th result as a picture (by default in "./out/")
component_index = 0
step = 5
save_plots(res, component_index, step)
```
## Generating Movies
The resulting images are stored (by default) in the `"./out/"` folder.
We can use ffmpeg to generate movies from these images.
```bash
ffmpeg \
	-y \
	-threads 16 \
	-r 24 \
	-f image2 \
	-pattern_type glob -i 'out/*.png' \
	-vcodec hevc \
	-pix_fmt yuv444p \
	-strict -2 \
	-tune animation \
	movie.mp4
```
Explaining the options:
```bash
# -y 									Overwrite output files
# -threads 16 							How many threads is ffmpeg going to use (this can safely be neglected)
# -r 24 								Framerate of the resulting movie
# -f image2								Format to use
# -pattern_type glob -i 'out/*.png'		Use the wildflag * to take many .png files as intput sources.
# -pix_fmt yuv444p						Pixel format which should work for png files saved with matplotlib.
# -strict -2							Specifies how strictly to follow standards.
# -tune animation						Tune the output for animations
```