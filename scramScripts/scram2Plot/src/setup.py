from setuptools import setup, find_packages

setup(
    name='scram2Plot',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
        'matplotlib>=3.5',
        'pandas>=1.4.0',
        'seaborn>=0.12.0',
        'numpy>=1.23.0',
    ],
    entry_points={
        'console_scripts': [
            'scram2Plot=scram2Plot:main',
        ],
        'scram2Plot': [
            'scram2Plot=scram2Plot:main',
        ]
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