"""
This is a setup script for installing the parser locally on python path with
all the required dependencies. Used mainly for local testing.
"""
from setuptools import setup, find_packages


#===============================================================================
def main():
    # Start package setup
    setup(
        name="nwchemparser",
        version="0.1",
        description="NoMaD parser implementation for NWChem.",
        author="Lauri Himanen",
        author_email="lauri.himanen@aalto.fi",
        license="GPL3",
        package_dir={'': 'parser/parser-nwchem'},
        packages=find_packages(),
        install_requires=[
            'pint',
            'numpy',
            'nomadcore',
        ],
    )

# Run main function by default
if __name__ == "__main__":
    main()
