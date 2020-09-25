"""
Build and install the package.
"""
from setuptools import find_packages, setup

NAME = 'pji'
FULLNAME = NAME
AUTHOR = "Matthias Steiner"
AUTHOR_EMAIL = 'matthias.steiner.@geo.tuwien.ac.at'
LICENSE = "BSD License"
URL = ""
DESCRIPTION = ""
LONG_DESCRIPTION = DESCRIPTION

VERSION = '1.0'

PACKAGES = find_packages()
SCRIPTS = []

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "Intended Audience :: Education",
    "Topic :: Scientific/Engineering",
    "Topic :: Software Development :: Libraries",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "License :: OSI Approved :: {}".format(LICENSE),
]
PLATFORMS = "Any"
INSTALL_REQUIRES = []  # see ../environment.yml

if __name__ == '__main__':
    setup(name=NAME, fullname=FULLNAME, description=DESCRIPTION,
          long_description=LONG_DESCRIPTION, version=VERSION, author=AUTHOR,
          author_email=AUTHOR_EMAIL, license=LICENSE, url=URL,
          platforms=PLATFORMS, scripts=SCRIPTS, packages=PACKAGES,
          classifiers=CLASSIFIERS, install_requires=INSTALL_REQUIRES)
