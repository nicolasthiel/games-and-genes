from setuptools import setup, find_packages

setup(
    name="games_and_genes",
    version="0.1.0",
    # This tells setuptools that package source code is in the 'src' directory
    package_dir={"": "src"},
    # This automatically finds packages inside 'src'
    packages=find_packages(where="src"),
    install_requires=[
        # List your dependencies here, e.g.,
        # "numpy",
        # "pandas",
    ],
)