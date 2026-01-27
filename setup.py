from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="virasign",
    version="1.0.0",
    author="Virasign Team",
    description="Virasign: Viral Read ASSIGNment from nanopore sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "virasign=virasign.virasign:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
