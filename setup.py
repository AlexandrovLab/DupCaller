from setuptools import find_packages, setup

setup(
    name="DupCaller",
    version="1.1.0",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/yuhecheng62/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=["src/DupCaller.py"],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        "biopython==1.85",
        "pysam==0.23.3",
        "numpy==2.3.4",
        "matplotlib==3.10.7",
        "scipy==1.16.2",
        "pandas==2.3.3",
        "h5py==3.15.0",
    ],
)
