# Create the Python runtime environment (legacy pyenv guide)

> **The project now uses uv.** See [`environment.md`](environment.md) for the
> current setup. This page is kept for reference if you prefer pyenv +
> pyenv-virtualenv; both honour the repo's `.python-version` (3.11.5).

## Background

The Python package dependencies are declared in
[`pyproject.toml`](../pyproject.toml) (with a committed `uv.lock`) and installed
with `uv sync`. See [`environment.md`](environment.md) for the current canonical
setup; this page covers the legacy pyenv + pip route.

This page goes through the Python environment setup steps in more detail.

**NOTE**: While you can create virtual environments using
`virtualenv` or `venv` in place of `pyenv`, be sure to put the environment
*outside* the `platelet-wcm/` directory.


## Prerequisites

* **Install** the software tools as described in [dev-tools](dev-tools.md), including
  * pyenv and pyenv-virtualenv
  * initialising pyenv in your shell profile
  * gcc or llvm
  * git
  * a programming editor such as PyCharm, Sublime Text, or Visual Studio Code
* **[Set up Git and GitHub](https://docs.github.com/en/get-started/quickstart/set-up-git)** including [Connecting to GitHub with SSH](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh).
* **Clone the repo** [platelet-wcm](https://github.com/stevehaigh/platelet-wcm) to a local directory like `~/dev/platelet-wcm/`. See [About remote repositories](https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories).


## Install native libraries

1. Use your package manager to install the needed libraries
   (most come from pyenv's requirements to install Python; see the
   [pyenv wiki](https://github.com/pyenv/pyenv/wiki) for the latest list).

   **On macOS**

   ```bash
   brew install openssl readline xz
   ```

   **On Ubuntu**

   ```bash
   sudo apt install -y libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev tk-dev
   ```


## Install Python

1. Use `pyenv`.

   ```bash
   pyenv install 3.11.5
   ```

   **Note:** If running the simulation gets an error related to "shared objects" and
   "needing to compile with -fPIC", use:

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.11.5
   ```


## Create the `platelet-wcm` python virtual environment

1. Create a pyenv virtualenv and select it in your project directory.

   ```bash
   cd ~/dev/platelet-wcm  # or wherever you cloned the repo
   pyenv virtualenv 3.11.5 platelet-wcm && pyenv local platelet-wcm
   ```

1. Upgrade this virtual environment's installers.

   ```bash
   pip install --upgrade pip setuptools wheel
   ```

1. Install the project dependencies (now declared in `pyproject.toml`; `uv sync`
   is the canonical install — see [`environment.md`](environment.md)). The repo is
   run from source (no install step), so for a bare-pip env install the pinned
   `[project.dependencies]` directly:

   ```bash
   pip install numpy==1.26.3 scipy==1.11.4 matplotlib==3.7.1 Unum==4.2.1 && pyenv rehash
   ```

1. **Required:** Add the following line to your shell profile and run it in your current shell.
   This gets more consistent results from OpenBLAS and improves performance significantly,
   especially when called from multiple processes.

   ```bash
   export OPENBLAS_NUM_THREADS=1
   ```

1. Set `PYTHONPATH` to the repo root. Consider adding an alias to your shell profile
   (see [dev-tools.md](dev-tools.md)):

   ```bash
   export PYTHONPATH=$PWD
   ```

1. Run the unit tests.

   ```bash
   python -m pytest models/platelet/tests/
   ```

1. Run a short simulation to verify the environment.

   ```bash
   python runscripts/manual/runPlateletSim.py out/test --length 60
   ```

   If either step fails with an error about "shared objects" or
   "needing to compile with -fPIC", reinstall Python with `--enable-shared`:

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.11.5 --force
   pyenv virtualenv-delete platelet-wcm
   # then repeat from the virtualenv creation step
   ```

1. If you're using PyCharm, select the project's Python interpreter so PyCharm
   understands the version of Python and its installed libraries. This enables
   code completion, visual debugging, click-through to library source code, etc.

   > PyCharm > 
   > Preferences > 
   > Project: platelet-wcm > 
   > Project Interpreter > 
   > gear ️ > 
   > Add... > 
   > Virtualenv Environment > 
   > Existing environment > 
   > Interpreter > 
   > [run `pyenv which python` in a shell to find the python location, something
   > like `/usr/local/var/pyenv/versions/platelet-wcm/python`, and paste that path
   > into the text box or navigate there].

