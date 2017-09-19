# GHER model postprocessing toolbox
A. Capet - Sept 2017

Main class is [G3D_class.py](G3D_class.py).
It contains the G3D class, and a serie of 

* variable definition : for variables derived from model outputs
* processing definition : integrate, average, along various dimensions, using masks, etc ... 
* utilities : eg. add derived variables to the model ouput files, ouptut new netcdf, plotting, etc ... 

Examples : 

* [EasyExample.py](Examples/EasyExample.py) (Annoted with explainations)
* [CentralVertProperties.py](Examples/CentralVertProperties.py)

User-specific variables ( in and out directories, configurations fiels , etc .. ) should be given in the file local.yml, for which a [template](local.yml.template) is provided 


For further documentation, refer to the [local wiki](https://github.com/acapet/GHER-POSTPROC/wiki)
