"""
Setup configuration for the Needleman-Wunsch Algorithm package.
"""

from setuptools import setup, find_packages

# Read the contents of README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="needleman-wunsch-aligner",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive Python implementation of the Needleman-Wunsch algorithm for global sequence alignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/needleman-wunsch-alignment",
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/needleman-wunsch-alignment/issues",
        "Documentation": "https://github.com/yourusername/needleman-wunsch-alignment/blob/main/README.md",
        "Source Code": "https://github.com/yourusername/needleman-wunsch-alignment",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    keywords="bioinformatics, sequence-alignment, needleman-wunsch, dynamic-programming, computational-biology",
    py_modules=["needleman_wunsch"],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.910",
        ],
        "examples": [
            "matplotlib>=3.3",
            "seaborn>=0.11",
        ],
    },
    entry_points={
        "console_scripts": [
            "needleman-wunsch=needleman_wunsch:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)