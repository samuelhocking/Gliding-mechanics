# Gliding mechanics

Author  : Sam Hocking
Updated : 5/10/2022

Final_project_Hocking.zip contains a few files:
- 1. Presentation.pdf : slide deck as presented in class
- 2. gliding.py : `python` file containing necessary classes and functions to implement the gliding simulations
- 3. gliding_notebook.ipynb : Jupyter notebook containing simulations of freefall, stable gliding, and dynamic gliding
- 4. YB49_cg_calc.nb : Mathematica notebook to approximate the x coordinate of the YB-49 center of gravity
- 5. atmosphere.csv : .CSV file containing data for the standard engineering atmosphere in SI units
- 6. naca653-018_coordinates.dat : coordinates in percentage of the chord for the NACA 653-018 airfoil
- 7. naca653018_polar_Re15k_clean.csv : polar for the NACA 653-018 airfoil at 1.5e4 Reynolds number
- 8. naca653018_polar_Re5mm_clean.csv : polar for the NACA 653-018 airfoil at 5e6 Reynolds number
- 9. naca653018_polar_Re10mm_clean.csv : polar for the NACA 653-018 airfoil at 1e7 Reynolds number
- 10. naca653018_polar_Re24mm_clean.csv : polar for the NACA 653-018 airfoil at 2.4e7 Reynolds number
## Installation

None required except dependent packages:
- numpy
- matplotlib
- pandas
- copy
- scipy

Place gliding.py and backup files 5-10 in the same directory as any notebook (including the attached notebooks) you wish to run with the contained classes and use the following import command:

```
from gliding import *
```
