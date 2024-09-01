Virtual Flask

Live corpus dashboard for publication: http://corpus.bmahjour.com/
Live analysis dashboard for publication: http://analysis.bmahjour.com/


For local use and development:

1. First, clone the repo:

git clone https://github.com/coleygroup/virtual_flask.git

2. Change directories:

cd virtual_flask

3. create a virtual environment:

python3 -m venv env

4. Activate it:

source env/bin/activate

5. Install the package

pip install -e .
(graphviz will have to be installed separtely, then pygraphviz, if network visualization code is desired)

6. Set bash script to executable

chmod +x bin/run

7. For an example notebook to shows how to use the flask and the filters, as well as the code for the figures in the manuscript, go to notebooks

8. To run the smarts toolkit:

bin/run corpus

and go to http://127.0.0.1:8000/ in your browser

9. To run the analysis software:

bin/run analysis

and go to http://127.0.0.1:8001/ in your browser

(you will need to have reduced pathway data in your local postgres server to use this dashboard)

Both apps can be run simultaneously. 

The templates can be found at shared/templates.



For advanced users, an HPC suite is contained in shared/hpc (QOL WIP)
For advanced users, npx is needed to work on the underlying javascript for dashboard visualizations