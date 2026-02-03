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
**reStructuredText (RST) style**, which is the standard for this project and 
provides seamless integration with Sphinx documentation.

Key guidelines:

1. **Function Docstrings** should describe the function purpose, parameters, return values, and exceptions.
2. **Class Docstrings** should describe the class purpose and provide an overview of its methods and attributes.
3. **Module Docstrings** should provide an overview of the module's functionality and main components.

RST Docstring Format
^^^^^^^^^^^^^^^^^^^^^

The RST style uses specific directives for documenting parameters, return values, and exceptions:

.. code-block:: python

    def example_function(param1: int, param2: str) -> bool:
        """
        Brief one-line summary of the function.

        More detailed description of the function's behavior, 
        implementation details, or usage notes can be added here
        in multiple paragraphs if needed.

        :param param1: Description of the first parameter
        :type param1: int
        :param param2: Description of the second parameter
        :type param2: str
        :return: Description of the return value
        :rtype: bool
        :raises NotImplementedError: Description of when this exception is raised
        :raises ValueError: Another exception that might be raised

        Example:
            >>> example_function(42, "hello")
            True
        
        .. note::
           Additional notes or warnings about the function.
        
        .. seealso::
           :func:`related_function`: Description of related functionality
        """
        raise NotImplementedError


    class ExampleClass:
        """
        Brief one-line summary of the class.

        Detailed description of the class purpose, main functionality,
        and usage patterns.

        :param init_param: Description of initialization parameter
        :type init_param: str
        :ivar attribute_name: Description of instance attribute
        :vartype attribute_name: int

        Example:
            >>> obj = ExampleClass("test")
            >>> obj.method()
            42
        
        .. note::
           Important information about class usage or behavior.
        """

        def __init__(self, init_param: str):
            """
            Initialize the ExampleClass instance.

            :param init_param: Parameter for initialization
            :type init_param: str
            """
            self.attribute_name: int = 0

        def method(self) -> int:
            """
            Brief description of the method.

            :return: Description of return value
            :rtype: int
            """
            return self.attribute_name


Common RST Directives
^^^^^^^^^^^^^^^^^^^^^

The following directives are commonly used in RST docstrings:

**Parameters and Return Values:**

- ``:param name:`` - Parameter description
- ``:type name:`` - Parameter type
- ``:return:`` - Return value description
- ``:rtype:`` - Return type

**Exceptions:**

- ``:raises ExceptionType:`` - Description of when exception is raised

**Class Attributes:**

- ``:ivar name:`` - Instance variable description
- ``:vartype name:`` - Instance variable type
- ``:cvar name:`` - Class variable description
- ``:vartype name:`` - Class variable type

**Cross-References:**

- ``:func:`function_name``` - Link to function
- ``:class:`ClassName``` - Link to class
- ``:mod:`module_name``` - Link to module
- ``:meth:`method_name``` - Link to method

**Admonitions:**

.. code-block:: restructuredtext

    .. note::
       Important information

    .. warning::
       Warning message

    .. seealso::
       Related references

    .. todo::
       Future improvements

**Examples:**

.. code-block:: restructuredtext

    Example:
        >>> code_example()
        result

Benefits of RST Style
^^^^^^^^^^^^^^^^^^^^^^

1. **Sphinx Integration**: RST docstrings are natively supported by Sphinx and render beautifully in the generated HTML documentation
2. **Type Documentation**: Explicit type directives improve code understanding
3. **Cross-Referencing**: Easy linking between different parts of documentation
4. **Flexibility**: Support for advanced formatting, admonitions, and examples
5. **Consistency**: Unified style across documentation and docstrings

.. tip::
   Use type hints in function signatures alongside RST docstrings for maximum clarity. 
   Sphinx can automatically extract type information from annotations, but explicit 
   ``:type:`` directives in docstrings improve readability.