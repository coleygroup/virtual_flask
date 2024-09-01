Virtual Flask



First, create a virtual environment:

python3 -m venv env

Activate it:

source env/bin/activate

Install the package

pip install -e .

To run the smarts toolkit:

bin/run corpus


To run the analysis software:

bin/run analysis


The templates can be found at shared/templates.


For an example notebook to shows how to use the flask and the filters, as well as the code for the figures in the manuscript, go to notebooks


For advanced users, an HPC suite is contained in shared/hpc (QOL WIP)