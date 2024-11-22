import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="WORMutils",
    version="0.0.3",
    author="Grayson Boyer",
    author_email="gmboyer@asu.edu",
    description="A package for common functions used in other WORM codes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    packages=['WORMutils'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['pandas', 'chemparse'],
    include_package_data=True,
    package_data={},
    zip_safe=False
)

