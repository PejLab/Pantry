import setuptools

setuptools.setup(
    name="Pantry",
    version="0.0.1",
    author="Daniel Munro",
    packages=["pantry"],
    install_requires=[
        "gtfparse",
        "numpy",
        "pandas",
    ]
)