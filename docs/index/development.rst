Development
===========


Project Structure Overview
--------------------------

.. toctree::
    :hidden:

The project is structured as shown below. Each folder serves a specific purpose 
which is described in the corresponding section where the main files and 
sub/folders are further described. 

.. code-block:: text

    QuditMeasurements/
    ├── .github/
    ├── .vscode/
    ├── .configs/
    ├── docs/
    ├── quaos/
    ├── scripts/
    └── tests/


.github
#######

This is the main folder for GitHub workflows. Its structured following the 
convention specified in the `GitHub Actions <https://docs.github.com/en/actions/writing-workflows/quickstart>`_
guide and its where the GitHub CI/CD pipelines are defined. 

.. code-block:: text

    QuditMeasurements/
    └── .github/
        └── workflows/


.vscode
#######

This is the main folder for Visual Studio Code configuration. 
It contains the definition of automated tasks that can be executed from the IDE; 
the suggested Visual Studio Code extension and their configuration. 

While, in general, it is not recommended to add this folder to version control, 
it is included here to speed up the configuration of the IDE and to provide a 
common starting point.

.. code-block:: text

    QuditMeasurements/
    └── .vscode/


.configs
########

In this folder are stored the configuration files for the project. This includes 
the .env file storing the project secrets (private keys, password, and so on). 

.. warning:: 
    The .env file is not versioned and should not be committed to the repository. 


.. code-block:: text

    QuditMeasurements/
    └── .configs/


docs
####

The docs folder contains source code necessary to build the project 
documentation. This includes both manually and automatically generated files 
and folders. The documentation 

.. code-block:: text

    QuditMeasurements/
    └── docs/       


scripts
#######

.. code-block:: text

    QuditMeasurements/
    └── scripts/
        ├── experiments/
        ├── notebooks/
        ├── personal/
        ├── python/
        └── shell/

quaos
#####

.. code-block:: text

    QuditMeasurements/
    └── quaos/ 


tests
#####

.. code-block:: text

    QuditMeasurements/
    └── tests/

    

Project Setup And Development Workflow
--------------------------------------

.. toctree::
    :hidden:

    development/tools
    development/extensions
    development/version_control
    development/testing
    development/documentation
    development/scripts


To set up the development environment for the project, follow the steps 
described in the guides below:

    1. `Required Development Tools <development/tools.html>`_ 
    2. `Visual Studio Code Extensions <development/extensions.html>`_  
    3. `Version Control <development/version_control.html>`_  
    4. `Testing <development/testing.html>`_
    5. `Project Documentation <development/documentation.html>`_
    6. `Scripts <development/scripts.html>`_


Coding Style
------------

The coding style adopted for this project is **PEP 8**. PEP 8 is the style guide 
for Python code, which helps to maintain readability and consistency across the 
codebase. While not mandatory, it is recommended to follow the this coding 
standard.

In `this page <https://www.python.org/dev/peps/pep-0008/>`_ 
is it possible to browse the full PEP 8 style guide. 

.. note::
   The suggested PEP 8-related extensions for Visual Studio Code (see 
   `Visual Studio Code Extensions <development/extensions.html>`_), automatically
   detect PEP 8 violations highlighting them in the editor as warnings
   and providing suggestions to fix them.
