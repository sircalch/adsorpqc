# Installation

AdsorpQC requires Python 3.9 or newer. Its dependencies are NumPy, SciPy, pandas, Matplotlib and Jinja2.

## From PyPI

```bash
pip install adsorpqc
```

Check the installation:

```bash
adsorpqc --version
```

## From source

```bash
git clone https://github.com/sircalch/adsorpqc.git
cd adsorpqc
pip install -e ".[dev]"
pytest
```

The `dev` extra installs the test and build tools.

## Building this documentation locally

```bash
pip install -e ".[docs]"
mkdocs serve
```

Then open <http://127.0.0.1:8000>.
