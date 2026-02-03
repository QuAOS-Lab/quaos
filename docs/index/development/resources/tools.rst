Required Development Tools
==========================

This guide outlines the key development tools used in our project. 
These tools help streamline our workflow, maintain code quality, 
and ensure smooth collaboration.


.. python:

Python
______

The codebase is primarily written in Python, leveraging its simplicity and 
versatility for development. The adopted version for this project is Python 3.11.9.

Windows
^^^^^^^

To install Python on Windows, it is recommended to 
download the installer form the `Python website <https://www.python.org/downloads/release/python-3119/>`_
(or by clicking `this link <https://www.python.org/downloads/release/python-3119/>`_) and run it.

Follow the installation wizard to complete the installation process. 
Once the installation is complete, it should be possible to run the following 
command from the command prompt:

.. code-block:: bash

    $ python --version


Linux
^^^^^

Python 3.11.9 can be installed on most Linux distributions using the package manager or by building from source.

**Ubuntu/Debian:**

.. code-block:: bash

    $ sudo apt update
    $ sudo apt install software-properties-common
    $ sudo add-apt-repository ppa:deadsnakes/ppa
    $ sudo apt update
    $ sudo apt install python3.11 python3.11-venv python3.11-dev

**Fedora/RHEL/CentOS:**

.. code-block:: bash

    $ sudo dnf install python3.11 python3.11-devel

**Arch Linux:**

.. code-block:: bash

    $ sudo pacman -S python

**Building from Source:**

If Python 3.11.9 is not available in your distribution's repository:

.. code-block:: bash

    $ wget https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz
    $ tar -xf Python-3.11.9.tgz
    $ cd Python-3.11.9
    $ ./configure --enable-optimizations
    $ make -j$(nproc)
    $ sudo make altinstall

Verify the installation:

.. code-block:: bash

    $ python3.11 --version


.. Git:

Git
___

Git is a distributed version control system used for tracking changes in source 
code. It allows developers to collaborate on a project, manage changes, and 
revert to previous versions. More Information about Git can be found 
`here <https://git-scm.com/>`_.

Windows
^^^^^^^

To install Git on Windows, it is recommended to 
download the last installer form the `Git website <https://git-scm.com/download/win>`_
and run it.

Once the installation is complete, it should be possible to run the following 
command from the command prompt:

.. code-block:: bash

    $ git --version


Linux
^^^^^

Git is typically available in most Linux distribution repositories and can be installed using the package manager.

**Ubuntu/Debian:**

.. code-block:: bash

    $ sudo apt update
    $ sudo apt install git

**Fedora/RHEL/CentOS:**

.. code-block:: bash

    $ sudo dnf install git

**Arch Linux:**

.. code-block:: bash

    $ sudo pacman -S git

Verify the installation:

.. code-block:: bash

    $ git --version

For the latest version, you can also install from the official Git PPA (Ubuntu/Debian):

.. code-block:: bash

    $ sudo add-apt-repository ppa:git-core/ppa
    $ sudo apt update
    $ sudo apt install git


.. Visual Studio Code:

Visual Studio Code
__________________

Visual Studio Code (VS Code) is a lightweight, powerful code editor that is 
highly customizable and widely used in the development community. With built-in 
support for debugging, Git integration, and a rich ecosystem of extensions, 
VS Code makes coding, testing, and version control more efficient.

.. note::
    While using Visual Studio Code is not mandatory, having a shared integrated 
    development environment (IDE) helps resolve common issues, improves teams
    efficiency and simplifies the development process.
    Moreover, all the other guides contained within this documentation assume 
    the use of VS Code as the IDE of choice.


Windows
^^^^^^^

To install Visual Studio Code on Windows, it is recommended to 
download the installer form the `Visual Studio Code website <https://code.visualstudio.com/download>`_
and run it.

Once the installation is complete, it should be possible to run the following 
command from the command prompt:

.. code-block:: bash

    $ code --version


Linux
^^^^^

Visual Studio Code can be installed on Linux through various methods depending on your distribution.

**Ubuntu/Debian (.deb package):**

Download and install the .deb package from the `VS Code website <https://code.visualstudio.com/download>`_:

.. code-block:: bash

    $ sudo apt install ./<file>.deb

Or install via the official Microsoft repository:

.. code-block:: bash

    $ wget -qO- https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor > packages.microsoft.gpg
    $ sudo install -D -o root -g root -m 644 packages.microsoft.gpg /etc/apt/keyrings/packages.microsoft.gpg
    $ sudo sh -c 'echo "deb [arch=amd64,arm64,armhf signed-by=/etc/apt/keyrings/packages.microsoft.gpg] https://packages.microsoft.com/repos/code stable main" > /etc/apt/sources.list.d/vscode.list'
    $ sudo apt update
    $ sudo apt install code

**Fedora/RHEL/CentOS (.rpm package):**

.. code-block:: bash

    $ sudo rpm --import https://packages.microsoft.com/keys/microsoft.asc
    $ sudo sh -c 'echo -e "[code]\nname=Visual Studio Code\nbaseurl=https://packages.microsoft.com/yumrepos/vscode\nenabled=1\ngpgcheck=1\ngpgkey=https://packages.microsoft.com/keys/microsoft.asc" > /etc/yum.repos.d/vscode.repo'
    $ sudo dnf check-update
    $ sudo dnf install code

**Arch Linux (AUR):**

.. code-block:: bash

    $ yay -S visual-studio-code-bin

Or using the snap package (works on most distributions):

.. code-block:: bash

    $ sudo snap install --classic code

Verify the installation:

.. code-block:: bash

    $ code --version

