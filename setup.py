from setuptools import setup, find_packages

setup(
    name='virtual_flask',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "sh",
        "wheel",
        "jinja2",
        "html5validator",
        "flask",
        "rdkit",
        "opencv-python",
        "pandas",
        "dill",
        "numpy<2",
        "psycopg2-binary",
        # "pygraphviz",
        "requests",
        "rdcanon @ git+https://github.com/coleygroup/rdcanon/"
    ],
    entry_points={
        'console_scripts': [
            'run-corpus=smarts_toolkit.webkit.__main__:main',
            'run-analysis=analysis.webkit.__main__:main',
            'run-campaign=campaign.webkit.__main__:main',
        ]
    },
)
