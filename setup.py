from setuptools import setup
from pyrpct.__init__ import __version__

setup(name='pyrpct',
    version=__version__,
    description='raac-pssm protein classification toolkit',
    url='https://github.com/KingoftheNight/pyrpct',
    author='liangyc',
    author_email='1694822092@qq.com',
    license='BSD 2-Clause',
    packages=['pyrpct'],
    install_requires=[
        'numpy',
        'matplotlib',
        'scikit-learn',
        'seaborn==-0.10.1',
        'pyecharts'
        ],
    entry_points={
        'console_scripts': [
        'pyrpct=pyrpct.__main__:rpct_main',
            ]
        },
    python_requires=">=3.6",
    include_package_data=True,
    zip_safe=True)
