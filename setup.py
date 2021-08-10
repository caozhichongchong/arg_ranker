from setuptools import setup

setup(
    name="arg_ranker",
    packages=['arg_ranker'],
    version="2.8.1",
    description="Ranking the risk of antibiotic resistance for genomes/metagenomes",
    author='An-Ni Zhang',
    author_email='anniz44@mit.edu',
    url='https://github.com/caozhichongchong/ARG_Ranker',
    keywords=['antibiotic resistance', 'risk', 'one health', 'clinical AMR', 'mobile AMR'],
    license='MIT',
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    package_dir={'arg_ranker': 'arg_ranker'},
    package_data={'arg_ranker': ['data/*','bin/*','example/*']},
    entry_points={'console_scripts': ['arg_ranker = arg_ranker.__main__:main']},
    install_requires=[
        'pandas',
        'argparse',
        'glob2',
        'statistics'
        ],
    classifiers=[
        "Programming Language :: Python :: 3"
    ]
)
