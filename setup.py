import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tigmint",
    version="1.2.8",
    author="Shaun Jackman",
    author_email="sjackman@gmail.com",
    description="Correct misassemblies using linked or long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bcgsc.github.io/tigmint/",
    license="GPLv3",
    python_requires=">=3",
    install_requires=[
        "intervaltree",
        "pybedtools",
        "pysam",
        "statistics",
        "numpy"],
    scripts=[
        "bin/tigmint",
        "bin/tigmint-arcs-tsv",
        "bin/tigmint-cut",
        "bin/tigmint-make",
        "bin/tigmint_molecule.py",
        "bin/tigmint_molecule_paf.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
