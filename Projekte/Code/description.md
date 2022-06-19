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
k1 = 10.0
k2 = 0.1
k3 = 4.93e-07
k4 = 80.0

def ODE(X, t):
    (A, B) = X
    return (k1 - k2*A + k3*A**2*B, k4 - k3*A**2*B)
```
To use the pde solver, we need to define initial values.
```python
initial_valus = np.random.uniform(low=200.0, high=900.0, size=(2, 10, 10))
```
Next, we initialize the solver and supply every information needed to obtain results.
```python
pde_solv = PDE_Solver(
	dx=40.0,
	dy=40.0,
	initial_values=initial_values,
	boundary_type="neumann",
	boundary_values=0.0,
	diffusion_constants=np.array([100.0, 5000.0]),
	kinetics=ODE
)
```
Now we can solve the PDE by specifying the times, we want to obtain.
```python
times = np.arange(0, 7200, 0.1)
res = pde_solv.solve_pde(times)
```
We can save plots by using the predefined method(s)
```python
# This will save the singular result at time step i of the results np.array
# The argument start_time can be neglected when calling individually
save_result_plot(i, res, index, min, max, start_time=None, output_folder=Path("./out/"))

# This will save every n_th result as a picture (by default in "./out/")
step = 500
index = 0
save_plots(res, index, step)
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