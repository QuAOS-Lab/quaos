Visual Studio Code Extensions
=============================

.. |python-badge| image:: https://img.shields.io/badge/Needed_for:-Python-blue
    :target: https://python.org
    :alt: For Python

.. |doc-badge| image:: https://img.shields.io/badge/Needed_for:-Docs-red
    :target: https://Sphinx-doc.org
    :alt: For Python

.. |github-badge| image:: https://img.shields.io/badge/Needed_for:-GitHub-gray
    :target: https://github.com
    :alt: For Python

To set up the Visual Studio Code configuration as described in the current workspace, follow these steps:

.. _install-extensions:

==========
Extensions
==========

To install a Visual Studio Code extension from the marketplace, follow these steps:  
   - Open the extensions view by clicking on the Extensions icon in the Activity 
     Bar on the side of the window or press ``Ctrl+Shift+X``.
   - Search for the extension by typing the name of the extension in the search 
     bar.
   - Click on the **Install** button to install the extension. 

The following extensions are strongly recommended for development: \

|python-badge|  
    - `Python extension pack <https://marketplace.visualstudio.com/items?itemName=donjayamanne.python-extension-pack>`_:    
      A collection of extensions for Python development.
    - `Flake8 <https://marketplace.visualstudio.com/items?itemName=ms-python.flake8>`_:  
      A linting tool for Python code.
    - `autopep8 <https://marketplace.visualstudio.com/items?itemName=ms-python.autopep8>`_:  
      Automatically formats Python code to conform to the PEP 8 style guide.
    - `Jupyter <https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter>`_:
      Provides Jupyter notebook support in Visual Studio Code.
    
|doc-badge| 
    - `reStructuredText <https://marketplace.visualstudio.com/items?itemName=lextudio.restructuredtext>`_:
      Provides reStructuredText language support.
    - `Esbonio <https://marketplace.visualstudio.com/items?itemName=swyddfa.esbonio>`_:
      Provides language server support for Sphinx documentation.

|github-badge|
    - `GitHub Actions <https://marketplace.visualstudio.com/items?itemName=GitHub.vscode-github-actions>`_:
      Provides support for GitHub Actions workflows.
    - `GitHub Copilot <https://marketplace.visualstudio.com/items?itemName=GitHub.copilot>`_:
      An AI-powered code completion tool.
    - `GitHub Copilot Chat <https://marketplace.visualstudio.com/items?itemName=GitHub.copilot-chat>`_:
      Provides chat-based interaction with GitHub Copilot.

The following extensions are optional:
    - `Material Icon Theme <https://marketplace.visualstudio.com/items?itemName=PKief.material-icon-theme>`_:
      A set of icons for Visual Studio Code.
    - `Open in Browser <https://marketplace.visualstudio.com/items?itemName=techer.open-in-browser>`_:
      Allows you to open files in your default web browser.
    - `Rainbow CSV <https://marketplace.visualstudio.com/items?itemName=mechatroner.rainbow-csv>`_:
      Provides syntax highlighting for CSV and TSV files.
    - `cpptools <https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools>`_:
      Provides C/C++ language support for Visual Studio Code.
    - `Task Runner <https://marketplace.visualstudio.com/items?itemName=SanaAjani.taskrunnercode>`_:
      Provides a task runner for Visual Studio Code.

1. **Configure Settings**:
   - Open the settings file located at [`.vscode/settings.json`](.vscode/settings.json) and ensure it contains the following configuration:

.. code-block:: json

    {
        "jupyter.notebookFileRoot": "${workspaceFolder}",
        
        "flake8.severity": 
        { 
            "C": "Information", 
            "E": "Warning", 
            "F": "Warning", 
            "R": "Hint", 
            "W": "Warning", 
            "I": "Information" 
        },
    }