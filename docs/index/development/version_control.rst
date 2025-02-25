Version Control
===============

This project uses Git for version control to track changes, collaborate 
effectively, and maintain a history of modifications. Below are guidelines for 
working with the repository.

In this section it is assumed that the readers is familiar with the 
basics of Git and how to use it. If this is not the case, please refer to the 
`Git documentation <https://git-scm.com/>`_ for more information.


Repository Structure
--------------------

The repository is organized in two main branches:    

- ``main``: contains the last stable version of the project.   

- ``dev``: contains the development version of the project with the most recent updates but potentially unstable.  


to ensures that the ``main`` branch remains stable and contains 
code that can be (eventually) distributed to third party collaborators 
(e.g. to researchers interested in the project that want to replicate some 
experiments) branch protection rules are put in place to prevent direct pushes. 
Only the ``dev`` branch can be merged with the ``main`` branch, and this 
happens when a new stable version of the project is released.

The ``dev`` branch is also protected by a branch protection rule that prevents 
direct pushes to the branch. The only way to update the ``dev`` branch is to 
merge a feature or bug fix back into it (more on this in `Git Workflow <#git_workflow>`_).  


.. git_workflow:

Git Workflow
------------

The general workflow for working on the project is as follows:

1. Clone the ``dev`` branch in your local machine.
2. Create a new feature-branch from local ``dev`` using a meaningful name that hints at the feature or bug fix you are working on.
3. Do some ✨magic✨ (i.e. implement the feature or fix the bug) on your feature-branch.
4. Push the changes to the remote repository creating (eventually) a corresponding remote feature-branch.
5. Open a merge request (MR) to merge the feature-branch back into ``dev``.
6. Wait for the approval of the MR by a code-owner.
7. (If required) Apply the suggested changes and go back to point 5.  
8. Once approved, merge the MR into ``dev``.

When working on a new feature or fixing a bug, you should create a new branch 
from ``dev`` and work on it. Once the feature or bug is complete, you can merge it 
back into ``dev``.


Cloning the Project
###################

To start working on the project, you need to clone the repository to your local 
machine. Use the following command:

.. code-block:: bash

    $ git clone -b dev <repository-url>

Replace `<repository-url>` with the actual URL of the Git repository. 
This will create a local copy of the ``dev`` branch of the project.


Create a Feature Branch
#######################

Once you have cloned the repository, you can start working on a new feature or 
bug fix. To create a new feature branch, use the following command:

.. code-block:: bash    

    $ git checkout -b <feature-branch-name>

Replace `<feature-branch-name>` with a meaningful name that hints at the feature 
or bug fix you are working on. This will create a new branch from the current 
state of the ``dev`` branch and "checkout" into it (i.e. setting it as the local 
working branch).

.. note:: 
    The checkout command is used to switch between branches. It is important to 
    keep track of which branch you are currently working on as the   


Apply Changes To The Feature Branch 
###################################

Once you have created a new feature branch, you can start working on the 
feature or bug fix. Apply your changes to the feature branch and commit them 
using the following command:

.. code-block:: bash

    $ git add .
    $ git commit -m "<commit-message>"

Make sure to replace `<commit-message>` with a clear and descriptive message 
that summarizes the changes you have made to the feature branch. This message 
will be included in the commit log when you push the changes to the remote 
repository.


Open a Merge Request
####################


