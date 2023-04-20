from setuptools import setup, find_packages

setup(
    name='adaptTrim',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'adaptTrim=adaptTrim:main',
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