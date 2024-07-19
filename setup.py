from setuptools import setup, find_packages

setup(
    name="prs-qc-toolkit",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pandas",
        "numpy",
        "pysam",
    ],
    author="Vamsee Pillalamarri",
    author_email="vpillal1@alumni.jh.edu",
    description="Quality control toolkit for Polygenic Risk Scores",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/vkp3/prs-qc-toolkit",
)
