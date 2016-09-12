This is the main repository of the [NOMAD](http://nomad-lab.eu) parser for
[NWChem](http://www.nwchem-sw.org/).

# Installation
This parser is a submodule of the nomad-lab-base repository. Developers within
the NoMaD project will automatically get a copy of this repository when they
download and install the base repository.

# Structure
The scala layer can access the parser functionality through the
scalainterface.py file, by calling the following command:

```python
    python scalainterface.py path/to/main/file
```

This scala interface is separated into it's own file to separate it from the
rest of the code. Some parsers will have the interface in the same file as the
parsing code, but I feel that this is a cleaner approach.

The parser is designed to support multiple versions of NWChem with a
[DRY](https://en.wikipedia.org/wiki/Don%27t_repeat_yourself) approach: The
initial parser class is based on NWChem 6.6, and other versions will be
subclassed from it. By sublassing, all the previous functionality will be
preserved, new functionality can be easily created, and old functionality
overridden only where necesssary.


# Standalone Mode
The parser is designed to be usable also outside the NoMaD project as a
separate python package. This standalone python-only mode is primarily for
people who want to easily access the parser without the need to setup the whole
"NOMAD Stack". It is also used when running custom unit tests found in the
folder *nwchem/test/unittests*. Here is an example of the call syntax:

```python
    from nwchemparser import NWChemParser
    import matplotlib.pyplot as mpl

    # 1. Initialize a parser by giving a path to the NWChem output file and a list of
    # default units
    path = "path/to/main.file"
    default_units = ["eV"]
    parser = NWChemParser(path, default_units=default_units)

    # 2. Parse
    results = parser.parse()

    # 3. Query the results with using the id's created specifically for NOMAD.
    scf_energies = results["energy_total_scf_iteration"]
    mpl.plot(scf_energies)
    mpl.show()
```

To install this standalone version, you need to clone the repositories
"python-common", "nomad-meta-info", and "parser-nwchem" into the same folder.
Then install the python-common according to the instructions found in the
README. After that, you can install this package by running either of the
following two commands depending on your python version:

```sh
python setup.py develop --user
python3 setup.py develop --user
```

# Tools and Methods
This section describes some of the guidelines that are used in the development
of this parser.

## Documentation
This parser tries to follow the [google style
guide](https://google.github.io/styleguide/pyguide.html?showone=Comments#Comments)
for documenting python code. Documenting makes it much easier to follow the
logic behind your parser.

## Testing
The parsers can become quite complicated and maintaining them without
systematic testing is impossible. There are general tests that are
performed automatically in the scala layer for all parsers. This is essential,
but can only test that the data is outputted in the correct format and
according to some general rules. These tests cannot verify that the contents
are correct.

In order to truly test the parser output, regression testing is needed. The
tests for this parser are located in
**/nwchem/parser/parser-nwchem/nwchemparser/regtest**. Tests provide one way to test
each parseable quantity and python has a very good [library for unit
testing](https://docs.python.org/2/library/unittest.html). When the parser
supports a new quantity it is quite fast to create unit tests for it. These
tests will validate the parsing, and also easily detect bugs that may rise when
the code is modified in the future.

## Profiling
The parsers have to be reasonably fast. For some codes there is already
significant amount of data in the NoMaD repository and the time taken to parse
it will depend on the performance of the parser. Also each time the parser
evolves after system deployment, the existing data may have to be reparsed at
least partially.

By profiling what functions take the most computational time and memory during
parsing you can identify the bottlenecks in the parser. There are already
existing profiling tools such as
[cProfile](https://docs.python.org/2/library/profile.html#module-cProfile)
which you can plug into your scripts very easily.
