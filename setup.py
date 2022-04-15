import setuptools

setuptools.setup(
    name="TURNAP",
    version="0.0.1",
    author="Daniel Munro",
    packages=["turnap"],
    install_requires=[
        "gtfparse",
        "numpy",
        "pandas",
    ]
)