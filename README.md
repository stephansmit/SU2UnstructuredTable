# Unstructured Thermodynamic Table

Code to generate unstructured thermodynamic VTK table using [CoolProp](http://www.coolprop.org/) and [Gmsh](http://gmsh.info/) 
<p align="center">
  <img src="https://github.com/stephansmit/SU2UnstructuredTable/blob/master/table.png" width=400>
</p>


## Requirements
Required Meshing software:
[Gmsh 4.4.0](http://gmsh.info/#Download)

Required Python modules:
```
python3 -m pip3 install --upgrade pip3 numpy
python3 -m pip3 install pandas pygmsh CoolProp meshio==3.0.0 matplotlib
```

## Use
Run locally
~~~
python3 make.table.py
~~~

Run using Singularity
~~~~
singularity pull shub://stephansmit/gmsh_containers
singularity exec gms_containers.sif python3 make_table.py
~~~~


