from setuptools import find_packages, setup

setup(
    name="DupCaller",
    version="1.2.1",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/AlexandrovLab/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=["src/DupCaller.py"],
    # src/ERROR (bundled fallback_latest.* error profiles) ships as the
    # data-only DupCaller_sub.ERROR package.
    package_dir={"": "src", "DupCaller_sub.ERROR": "src/ERROR"},
    packages=find_packages(where="src") + ["DupCaller_sub.ERROR"],
    package_data={"DupCaller_sub.ERROR": ["fallback_latest.*.txt"]},
    install_requires=[
        "biopython==1.85",
        "pysam==0.23.3",
        "numpy==2.3.4",
        "matplotlib==3.10.7",
        "scipy==1.16.2",
        "pandas==2.3.3",
        "h5py==3.15.0",
        "sigProfilerPlotting==1.4.3",
    ],
    extras_require={
        "test": ["pytest"],
    },
)
