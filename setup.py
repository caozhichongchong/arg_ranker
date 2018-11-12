from setuptools import setup

setup(
    name="arg_ranker",
    packages=['arg_ranker'],
    version="1,0",
    description="Ranking the risk of antibiotic resistance for metagenomes",
    author='Anni Zhang',
    author_email='anniz44@mit.edu',
    url='https://github.com/caozhichongchong/ARG_Ranker',
    keywords=['antibiotic resistance', 'risk', 'one health', 'clinical AMR', 'mobile AMR'],
    install_requires=['python>=3.0'],
    #include_package_data=True,
    long_description=open('README.md').read(),
    packages=find_packages(exclude=['argparse', 'os', 'csv']),
    package_data={  # Optional
        'ARG_rank': ['arg_ranker/ARG_rank.txt'],
    },
    #entry_points={'console_scripts': ['pyinstrument = pyinstrument.__main__:main']},
    #zip_safe=False,
    #setup_requires=['pytest-runner'],
    #tests_require=['pytest'],
    classifiers=[
        'Development Status :: 1 - Aplha',
        #'Environment :: Console',
        #'Environment :: Web Environment',
        'Intended Audience :: Bioinformatics and Researchers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Antibiotic resistance :: risk ranking',
        'Topic :: Metagenomes :: Antibiotic resistance',
    ]
)
