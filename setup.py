from setuptools import setup, find_packages


def main():
    setup(
        name="nwchemparser",
        version="0.1",
        description="NOMAD parser implementation for NWChem.",
        author="Lauri Himanen",
        author_email="lauri.himanen@aalto.fi",
        license="GPL3",
        package_dir={'': 'parser/parser-nwchem'},
        packages=find_packages(),
        install_requires=[
            'nomadcore',
        ],
    )

if __name__ == "__main__":
    main()
