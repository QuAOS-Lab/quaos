.. _version_control:

=========================
Version Control & Gitflow
=========================

This project uses Git for version control following a structured branching model to ensure stability, traceability, and seamless collaboration. This document outlines the branching strategy, versioning scheme, and workflow guidelines.

.. note::
   This guide assumes familiarity with Git basics. If you're new to Git, please refer to the `official Git documentation <https://git-scm.com/doc>`_ before proceeding.


Branching Model
===============

The repository follows a protected dual-branch model with strict merge policies enforced through GitHub branch protection rules.

Core Branches
-------------

**dev Branch**
   The ``dev`` branch is the **main development branch** and serves as the integration point for all features and bug fixes. 
   
   - Acts as the primary branch for active development
   - All feature branches merge into ``dev`` via pull requests
   - Every commit to ``dev`` is automatically tagged with semantic versioning
   - Protected branch: direct pushes are forbidden
   - Commits can only be added through approved pull requests

**main Branch**
   The ``main`` branch is the **production-ready branch** containing only stable releases.
   
   - Represents the latest stable release (follows ``Major.Minor.Patch`` versioning)
   - Contains thoroughly tested and approved code suitable for distribution
   - Protected branch: direct pushes are forbidden
   - Updates only through dedicated release branches or pull requests
   - Intended for external collaborators and researchers requiring stable code

.. important::
   The ``main`` and ``dev`` branches **cannot be merged directly into each other**. Specific update branches must be created to synchronize changes between these branches when necessary.


Versioning Scheme
=================

The project follows `Semantic Versioning 2.0.0 <https://semver.org/>`_ (SemVer) with automated tagging on the ``dev`` branch.

Version Format
--------------

Versions follow the ``MAJOR.MINOR.PATCH`` format:

- **MAJOR**: Incremented for incompatible API changes
- **MINOR**: Incremented for backward-compatible functionality additions
- **PATCH**: Incremented for backward-compatible bug fixes

**Development Versions:**

Commits to ``dev`` are automatically tagged with development version identifiers that include distance from the last stable release:

.. code-block:: text

   v0.2.1.12

Where:
   - ``v0.2.1``: Last stable release version
   - ``12``: 12 commits since last stable release

**Release Versions:**

Releases on ``main`` use clean version tags without development suffixes:

.. code-block:: text

   v1.0.0
   v1.2.3

Automated Tagging
-----------------

The repository uses GitHub Actions to automatically tag commits:

- **dev Branch**: Every commit triggers the ``SemVer Tagging`` workflow (`.github/workflows/tag_commit.yml`) which determines and applies the appropriate development version tag
- **main Branch**: Release versions are created through GitHub releases and tagged manually following the SemVer format

.. note::
   Automated tagging ensures complete version traceability and enables proper version resolution during package builds.


Branch Protection Rules
=======================

Both ``main`` and ``dev`` branches have the following protections enabled:

1. **No Direct Pushes**: All changes must be submitted via pull requests
2. **Required Reviews**: Pull requests require approval from code owners before merging
3. **Status Checks**: Automated CI/CD workflows must pass before merging
4. **Linear History**: Enforces a clean commit history through squash or rebase merging

These protections ensure code quality, prevent accidental overwrites, and maintain a reviewable history of all changes.


Git Workflow
============

The recommended workflow follows these steps to maintain branch integrity and enable collaborative development.

1. Clone the Development Branch
--------------------------------

Start by cloning the ``dev`` branch to your local machine:

.. code-block:: bash

   $ git clone -b dev <repository-url>
   $ cd <repository-name>

Replace ``<repository-url>`` with the actual repository URL (available on the GitHub repository page).

2. Create a Feature Branch
---------------------------

Create a new feature branch from ``dev`` using a descriptive name:

.. code-block:: bash

   $ git checkout -b feature/your-feature-name

**Branch Naming Conventions:**

- ``feature/``: New features or enhancements (e.g., ``feature/pauli-symmetry-detection``)
- ``bugfix/``: Bug fixes (e.g., ``bugfix/circuit-measurement-error``)
- ``hotfix/``: Critical fixes requiring immediate attention
- ``refactor/``: Code restructuring without changing functionality
- ``docs/``: Documentation updates
- ``test/``: Adding or updating tests

.. tip::
   Use descriptive, lowercase names with hyphens separating words. Include relevant issue numbers when applicable (e.g., ``feature/issue-42-optimize-hamiltonian``).

3. Implement Changes
--------------------

Make your changes on the feature branch, committing incrementally with clear messages:

.. code-block:: bash

   $ git add <changed-files>
   $ git commit -m "Brief description of changes"

**Commit Message Guidelines:**

- Use imperative mood ("Add feature" not "Added feature")
- Keep the first line under 50 characters
- Provide detailed explanation in the body if needed
- Reference related issues with ``#issue-number``

Example:

.. code-block:: bash

   $ git commit -m "Add Pauli symmetry detection algorithm
   
   Implements symmetry detection for stabilizer groups using
   symplectic representation. Resolves #42."

4. Push Feature Branch
----------------------

Push your feature branch to the remote repository:

.. code-block:: bash

   $ git push -u origin feature/your-feature-name

This creates a corresponding remote branch and sets up tracking.

5. Open a Pull Request
----------------------

On GitHub, open a pull request (PR) to merge your feature branch into ``dev``:

1. Navigate to the repository on GitHub
2. Click "Pull requests" → "New pull request"
3. Set base branch to ``dev`` and compare branch to your feature branch
4. Fill in the PR template with:
   
   - Clear description of changes
   - Related issue numbers
   - Testing performed
   - Breaking changes (if any)

5. Request review from appropriate code owners

.. note::
   The PR will automatically trigger CI/CD workflows defined in ``.github/workflows/`` that perform code quality checks, run tests, and verify build integrity.

6. Address Review Feedback
--------------------------

If reviewers request changes:

.. code-block:: bash

   $ git add <modified-files>
   $ git commit -m "Address review feedback: <description>"
   $ git push

The PR will automatically update with your new commits.

7. Merge the Pull Request
--------------------------

Once approved and all checks pass, merge the PR:

- Use "Squash and merge" for feature branches to keep ``dev`` history clean
- Use "Rebase and merge" if commit history on the feature branch is well-structured
- Delete the feature branch after merging to keep the repository tidy

After merging, the ``SemVer Tagging`` workflow automatically creates a new version tag on ``dev``.


Release Workflow
================

Creating a stable release on ``main`` requires careful coordination to ensure quality and proper versioning.

Creating a Release
------------------

1. **Prepare Release Branch**

   Create a release preparation branch from ``dev``:

   .. code-block:: bash

      $ git checkout dev
      $ git pull origin dev
      $ git checkout -b release/v1.2.0

2. **Final Testing and Documentation**

   - Update version numbers in relevant files
   - Update CHANGELOG with release notes
   - Run comprehensive test suites
   - Build and verify documentation

3. **Create Release PR to main**

   Open a pull request to merge the release branch into ``main``:

   - Ensure all tests pass
   - Get approval from project maintainers
   - Verify version tags are correct

4. **Create GitHub Release**

   After merging to ``main``, create a GitHub release:

   - Navigate to "Releases" → "Draft a new release"
   - Choose the appropriate tag (e.g., ``v1.2.0``)
   - Write release notes highlighting changes
   - Mark as pre-release if applicable
   - Publish release

   This triggers the ``Build and Publish Wheels`` workflow, which builds Python packages and publishes to PyPI.

5. **Synchronize Changes Back to dev**

   After releasing, create a branch to merge changes from ``main`` back to ``dev``:

   .. code-block:: bash

      $ git checkout main
      $ git pull origin main
      $ git checkout -b sync/main-to-dev
      $ git merge main
      $ git push -u origin sync/main-to-dev

   Open a PR to merge ``sync/main-to-dev`` into ``dev``.


Hotfix Workflow
---------------

For critical bugs in production (``main``):

1. Create hotfix branch from ``main``:

   .. code-block:: bash

      $ git checkout main
      $ git checkout -b hotfix/critical-bug-fix

2. Implement fix and create PR to ``main``

3. After merging to ``main``, create release with patch version bump

4. Sync hotfix back to ``dev`` using sync branch


Special Branches
================

Beyond the core ``dev`` and ``main`` branches, several special-purpose branches may exist:

Release Branches
----------------

**Format**: ``release/vX.Y.Z``

Created when preparing a new release. Allows final testing and version bumping without blocking ongoing development on ``dev``.

Sync Branches
-------------

**Format**: ``sync/source-to-target``

Used to synchronize changes between ``main`` and ``dev`` since direct merges are forbidden. Examples:

- ``sync/main-to-dev``: Brings release changes and hotfixes from ``main`` to ``dev``
- ``sync/dev-to-main``: Prepares new release from ``dev`` to ``main``

Update Branches
---------------

**Format**: ``update/dependency-name`` or ``update/area``

Used for dependency updates, security patches, or systematic refactoring across the codebase.


Best Practices
==============

1. **Keep Branches Short-Lived**
   
   Merge feature branches quickly to avoid conflicts and integration issues.

2. **Sync Regularly**

   Pull latest changes from ``dev`` frequently:

   .. code-block:: bash

      $ git checkout dev
      $ git pull origin dev
      $ git checkout feature/your-feature
      $ git merge dev

3. **Write Meaningful Commits**

   Each commit should represent a logical unit of change with a clear message.

4. **Test Before Pushing**

   Run local tests before pushing to avoid breaking CI/CD pipelines.

5. **Review Your Own Code First**

   Review your diff before opening a PR to catch obvious issues.

6. **Keep PRs Focused**

   One PR should address one feature or bug. Split large changes into multiple PRs.

7. **Update Documentation**

   Include documentation updates in the same PR as code changes.


Continuous Integration
======================

The repository uses GitHub Actions for automated workflows:

``check_source.yml``
   Validates code quality, runs linters, and performs static analysis

``review_code.yml``
   Executes test suite and generates coverage reports

``tag_commit.yml``
   Automatically tags commits on ``dev`` with SemVer development versions

``publish_repo.yml``
   Builds Python wheels and publishes to PyPI on release creation

All workflows must pass before a PR can be merged.


Troubleshooting
===============

Cannot Push to dev or main
---------------------------

Both branches are protected. Create a feature branch and submit a PR instead.

Merge Conflict on PR
--------------------

Update your branch with latest ``dev``:

.. code-block:: bash

   $ git checkout feature/your-feature
   $ git fetch origin
   $ git merge origin/dev
   # Resolve conflicts
   $ git add <resolved-files>
   $ git commit
   $ git push

PR Checks Failing
-----------------

Review the GitHub Actions logs to identify failures. Fix issues locally, commit, and push to update the PR.

Need to Update main from dev
-----------------------------

Do not merge directly. Create a ``release/vX.Y.Z`` branch, open a PR to ``main``, and follow the release workflow.


See Also
========

- `Git Documentation <https://git-scm.com/doc>`_
- `Semantic Versioning <https://semver.org/>`_
- `GitHub Flow Guide <https://guides.github.com/introduction/flow/>`_
- :ref:`development_tasks`: Automated development tasks



