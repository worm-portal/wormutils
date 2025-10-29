import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wormutils",
    version="0.2.0",
    author="Grayson Boyer",
    author_email="gmboyer@asu.edu",
    description="A package for common functions used in other WORM codes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    packages=['wormutils'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=['pandas', 'chemparse', 'requests', 'matplotlib', 'numpy', 'importlib_resources'],
    include_package_data=True,
    package_data={},
    zip_safe=False
)

