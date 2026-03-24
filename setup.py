from setuptools import find_packages, setup
from hapi.hapi import HAPI_VERSION

setup(
    name="hitran-api",
    version=HAPI_VERSION,
    packages=find_packages(include=["hapi", "hapi.*", "research", "research.*"]),
    license="MIT",
)
