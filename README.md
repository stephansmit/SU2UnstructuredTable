
# Unstructured Thermodynamic Table
<img align="right" width="250" src="https://github.com/stephansmit/SU2UnstructuredTable/blob/master/table.png">

Code to generate unstructured thermodynamic VTK table using [CoolProp](http://www.coolprop.org/) and [Gmsh](http://gmsh.info/) 

## Requirements
Required Meshing software:
[Gmsh 4.4.0](http://gmsh.info/#Download)

Required Python modules:
```
pip3 install --upgrade pip3 numpy
pip3 install pandas pygmsh CoolProp meshio==3.0.0 matplotlib
```

Or use Singularity container, see below.
## Use
Run locally
~~~
python3 make.table.py
~~~

Run using [Singularity](https://sylabs.io/singularity/) 
~~~~
singularity pull shub://stephansmit/gmsh_containers
singularity exec gms_containers.sif python3 make_table.py
~~~~


