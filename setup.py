import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tigmint",
    version="1.1.2",
    author="Shaun Jackman",
    author_email="sjackman@gmail.com",
    description="Correct misassemblies using linked reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bcgsc.github.io/tigmint/",
    license="GPLv3",
    python_requires=">=3",
    install_requires=[
        "intervaltree",
        "pybedtools",
        "pysam",
        "statistics"],
    scripts=[
        "bin/tigmint",
        "bin/tigmint-arcs-tsv",
        "bin/tigmint-cut",
        "bin/tigmint-make",
        "bin/tigmint-molecule"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
