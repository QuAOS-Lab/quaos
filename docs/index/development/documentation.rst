Project Documentation
=====================

Our project uses an automated documentation generation tool. 
To build the documentation, run these commands in the terminal:

.. code-block:: bash

    $ cd docs
    $ make html


Docstrings
----------

To ensure that the documentation is comprehensive and up-to-date, all functions 
and classes must include docstrings. These docstrings should follow the 
`Google Python Style Guide <https://google.github.io/styleguide/pyguide.html>`_.
In particular:

1. **Function Docstrings** should describes the function purpose, parameters, and return values.
2. **Class Docstrings** should describes the class purpose and provides an overview of its methods and attributes.


.. code-block:: python

    def example_function(param1: int, param2: str) -> bool:
        """
        This function demonstrates how to write a docstring.

        Args:
            param1 (int): The first parameter.
            param2 (str): The second parameter.

        Returns:
            bool: The return value. True for success, False otherwise.
        
        Raises:
            NotImplementedError: If the function is not implemented.
        """
        raise NotImplementedError