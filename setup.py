from setuptools import setup

setup(
    name="traits_finder",
    packages=['traits_finder'],
    version="1.5",
    description="search and summarize traits in genomes and metagenomes",
    author='Anni Zhang',
    author_email='anniz44@mit.edu',
    url='https://github.com/caozhichongchong/traits_finder',
    keywords=['metagenomes', 'genomes', 'function', 'traits'],
    license='MIT',
    install_requires=[
    'biopython',
    'argparse',
    'glob2',
    'statistics'
    ],
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    package_dir={'traits_finder': 'traits_finder'},
    package_data={'traits_finder': ['scripts/*','database/*','example/*','*.py']},
    entry_points={'console_scripts': ['traits_finder = traits_finder.__main__:main']},
    classifiers=[
        #'Development Status :: 1 - Alpha',
        #'Intended Audience :: Bioinformatics and Researchers',
        #'License :: MIT',
        #'Operating System :: MacOS',
        #'Operating System :: Microsoft :: Windows',
        #'Operating System :: LINUX',
        'Programming Language :: Python :: 3',
        #'Topic :: Antibiotic resistance :: risk ranking',
        #'Topic :: Metagenomes :: Antibiotic resistance',
    ]
)
