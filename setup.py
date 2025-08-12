from setuptools import *

setup(
    name="DupCaller",
    version="1.0.5",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/yuhecheng62/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=["src/DupCaller.py"],
    package_dir={"": "src"},
    install_requires=[
        "biopython>=1.78",
        "pysam>=0.19.0",
        "numpy>=1.21.5",
        "matplotlib>=3.8.0",
        "scipy>=1.11.3",
        "pandas>=1.5.3",
        "h5py>=3.10.0",
    ],
)
