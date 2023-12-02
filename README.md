# Y2 Computing Project

## Overview

This repository, `Y2Computing-Project`, is dedicated to housing one of the two Y2 Imperial College Physics computing projects focused on Python-based Object-Oriented development. 

## Repository Contents

- `.idea/`: PyCharm project files.
- `.gitignore`: Lists files and directories to be ignored by git.
- `plot_generator.py`: A Python script for generating plots for the report.
- `raytracer.py`: Where all of the physics-related python is housed. It has a number of classes that are fairly well documented so its use should be relatively easy.
- `testfile.py`: Python script for testing.

## Usage

The primary files of interest are `plot_generator.py`, `raytracer.py`, and `testfile.py`. These can be run in any Python environment.
All packages are included in the standard python install.

raytracer.py will raise warnings in some occasions (when ray is terminated, usually). This is particularly useful when
debugging and testing, so even though it might clutter the terminal somewhat, should someone else develop this code
further, they would appreciate this feature.

plot_generator.py might seem like an very poorly written script, however it is like this by design. The code is
meant to able to run by cells, by simply running the first cell with the import statements, then the rest individually. 
Furthermore, the explicit code repetition makes it especially easy to see what is going on for
testing and debugging.

In plot_generator.py all plots run quite quickly, except for the last one, since it has to propagate a very large number
of rays. If you just want to see the functionality then it might be worth changing the following parameters so less
datapoints are obtained, and the code runs faster:
    In lines 106, 107, 114, 115, 122 and 123 the parameters with value 0.3 may be changed to a larger value (e.g. 0.5)

## License

This project is under the MIT License - see the [LICENSE.md](LICENSE.md) file for details, including the requirement to cite this work in any derivative uses. If you are a current Y2 Imperial Physics student, you shouldn't use this code, probably.
