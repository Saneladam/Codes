# MICROMAMBA
Usa micromamba para:
*OpenMC, ML, CUDA, HPC, librerías difíciles, stacks científicos.

``` bash
# to create
micromamba create -n openmc python=3.12
# to activate
micromamba activate openmc
# to set
micromamba install -c conda-forge openmc
# to deactivate
micromamba deactivate
```

# PYTHON-VENV
Usa venv para:
*scripts normales, proyectos pequeños, web, automatización, tooling.

``` bash
# to create
python -m venv .venv
# to activate
source .venv/bin/activate
# to set
pip install numpy matplotlib
# to deactivate
deactivate
```
