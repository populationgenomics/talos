"""
automated installation instructions
"""


from setuptools import find_packages, setup


with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list

    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [
            line.strip()
            for line in filehandler
            if line.strip() and not line.startswith('#')
        ]


setup(
    name='automated-interpretation-pipeline',
    description='CPG Variant Prioritisation',
    long_description=readme,
    version='2.0.3',
    author='Matthew Welland, CPG',
    author_email=(
        'matthew.welland@populationgenomics.org.au, '
        'cas.simons@populationgenomics.org.au'
    ),
    package_data={'reanalysis': ['templates/*.jinja']},
    url='https://github.com/populationgenomics/automated-interpretation-pipeline',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=read_reqs('requirements.txt'),
    extras_require={
        'test': read_reqs('requirements-dev.txt'),
    },
)
