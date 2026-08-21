from setuptools import setup, find_packages

setup(
    name='scram2Plot',
    version='0.0.3',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.9',
    install_requires=[
        'matplotlib>=3.5',
        'pandas>=1.4.0',
        'numpy>=1.23.0',
    ],
    extras_require={
        # secondary-structure plots (profilePlot --sec_structure, utils MFE)
        'structure': ['biopython>=1.79', 'viennarna>=2.4.18'],
        # sequence logos (utils.ScramCSV.generate_sequence_logo)
        'logo': ['logomaker'],
        # parallel dataframe engine for very large CSVs (utils.ScramCSV(engine="modin"))
        'modin': ['modin', 'grpcio'],
    },
    entry_points={
        'console_scripts': [
            'profilePlot=scram2Plot.profilePlot:main',
            'binPlot=scram2Plot.binPlot:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
