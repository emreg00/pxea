from setuptools import setup

setup(name='pxea',
    version='1.0',
    description='Proximal pathway enrichment analysis',
    long_description='Methods for running ProXimal pathway Enrichment Analysis (PXEA). See project web page for more info.',
    url='http://github.com/emreg00/pxea',
    author='Emre Guney',
    author_email='emreguney@gmail.com',
    license='',
    packages=['pxea'],
    install_requires=[
	          'numpy',
	          'scipy',
	          'networkx',
		    ],
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe=False)

