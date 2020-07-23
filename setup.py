import setuptools

setuptools.setup(
    name="pymolsim",
    version="0.0.1",
    author="Sarath Menon, Grisell Diaz Leines, Jutta Rogal",
    author_email="sarath.menon@rub.de",
    description="Simple toy md code",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages('src'),
    package_dir={'':'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #install_requires=['numpy', 'joblib'],
    #scripts=['bin/pathsampling', 'bin/pathsampling_kernel'],
)
