from setuptools import setup, find_packages

setup(
    name='scram2Plot',
    version='0.0.2',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
        'matplotlib>=3.5',
        'pandas>=1.4.0',
        'numpy>=1.23.0',
        'biopython>=1.79',
        'viennarna>=2.4.18',
        'logomaker',
        'modin',
        'ray>=2.7.0',
        'grpcio'

    ],
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
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)