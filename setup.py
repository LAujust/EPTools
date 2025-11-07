from setuptools import setup, find_packages

setup(
    name="EPTools",                     # package name
    version="0.0.1",                      # initial version
    author="R. -D. Liang",
    author_email="liangrd@bao.ac.cn",
    description="An easy package/code for Einstein Probe (or other, e.g. XRT) data fetching and analysis.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/myproject",  # optional
    packages=find_packages(),              # automatically finds all packages
    install_requires=[
        "numpy>=2.3",
        "matplotlib>=3.5",
        "pandas>=2.2.2"
    ],
    python_requires=">=3.8",               # optional
)