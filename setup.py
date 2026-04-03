from setuptools import setup, find_packages

setup(
    name="geneXref",
    version="0.1.0",
    packages=find_packages(),
    package_data={"geneXref": ["data/*.tsv"]},
    install_requires=["pandas"],
    extras_require={
        "dev": ["pytest"],
    },
)
