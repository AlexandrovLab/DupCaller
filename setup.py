from setuptools import find_packages, setup

setup(
    name="DupCaller",
    version="1.2.3",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/AlexandrovLab/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=["src/DupCaller.py"],
    # src/ERROR (bundled fallback_latest.* error profiles) ships as the
    # data-only DupCaller_sub.ERROR package.
    package_dir={"": "src", "DupCaller_sub.ERROR": "src/ERROR"},
    packages=find_packages(where="src") + ["DupCaller_sub.ERROR"],
    # DupCaller_sub.PERF is a vendored copy of PERF (https://github.com/rkmlab/perf,
    # MIT license, see src/DupCaller_sub/PERF/LICENSE); its lib/ assets are
    # needed for PERF's HTML report (-a).
    package_data={
        "DupCaller_sub.ERROR": ["fallback_latest.*.txt"],
        "DupCaller_sub.PERF": [
            "LICENSE",
            "README.md",
            "all_repeats_1-6nt.txt",
            "lib/*.html",
            "lib/src/*.js",
            "lib/styles/*.css",
        ],
    },
    entry_points={"console_scripts": ["PERF=DupCaller_sub.PERF.core:main"]},
    install_requires=[
        "biopython==1.85",
        "pysam==0.23.3",
        "numpy==2.3.4",
        "matplotlib==3.10.7",
        "scipy==1.16.2",
        "pandas==2.3.3",
        "h5py==3.15.0",
        "sigProfilerPlotting==1.4.3",
        "tqdm>=4",
    ],
    extras_require={
        "test": ["pytest"],
    },
)
