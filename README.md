# Virtual Flask: Local Development Guide

## Overview

The Virtual Flask is a cheminformatic application designed to generate reaction networks given a set of initial molecules and a corpus of templates. A suite of supporting tools facilitate the analysis and visualization of generated chemical reaction data. This guide provides detailed instructions for setting up and running the Virtual Flask locally for development purposes. Both the Corpus and Analysis dashboards can be hosted locally and run simultaneously.

### Live Corpus Dashboards

- [Live Corpus Dashboard](http://corpus.bmahjour.com/)
- [Live Analysis Dashboard](http://analysis.bmahjour.com/)

---

## Local Setup Instructions

Follow the steps below to set up and run Virtual Flask locally for development purposes.

### **1. Clone the Repository**

First, clone the Virtual Flask repository to your local machine:

```sh
git clone https://github.com/coleygroup/virtual_flask.git
```

### **2. Navigate to the Project Directory**

Change directories to enter the project folder:

```sh
cd virtual_flask
```

### **3. Create a Virtual Environment**

Create a Python virtual environment to manage dependencies:

```sh
python3 -m venv env
```

### **4. Activate the Virtual Environment**

Activate the virtual environment:

- On Linux/macOS:
  ```sh
  source env/bin/activate
  ```
- On Windows (PowerShell):
  ```sh
  .\env\Scripts\Activate.ps1
  ```

### **5. Install Dependencies**

Install the package in editable mode:

```sh
pip install -e .
```

> **Note**: To enable network visualizations, you will need to install `graphviz` separately and then install `pygraphviz`.

### **6. Set Bash Script to Executable**

Make the script that runs the app executable:

```sh
chmod +x bin/run
```

### **7. Explore Example Notebooks**

To learn how to use the Flask application and its filtering features, as well as review the code used for figures in the related manuscript, explore the provided Jupyter notebooks:

- Navigate to the `notebooks` directory.

### **8. Run the Corpus Dashboard**

To start the Corpus Dashboard locally:

```sh
bin/run corpus
```

Then, open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) in your browser to view the app.

### **9. Run the Analysis Dashboard**

To start the Analysis Dashboard:

```sh
bin/run analysis
```

Open [http://127.0.0.1:8001/](http://127.0.0.1:8001/) in your browser to view the analysis app.

> **Note**: You will need reduced pathway data stored in your local PostgreSQL server to fully utilize this dashboard.

### **10. Running Both Apps Simultaneously**

Both the Corpus and Analysis dashboards can be run at the same time without any conflicts.

---

## Additional Resources

### **Templates**

- HTML templates for the dashboards are located in the `shared/templates` directory.

### **Advanced User Features**

- **HPC Suite**: An HPC suite is included in `shared/hpc`. This suite allows for the execution of high-performance computational workflows. To use it:

  1. Navigate to the `shared/hpc` directory:

     ```sh
     cd shared/hpc
     ```

  2. Customize the available scripts to fit your workload requirements. The scripts are designed to help automate and streamline batch processing on high-performance computing clusters.

  3. Run the desired script using a compatible job scheduler, such as `SLURM` or `PBS`.

  > **Note**: The HPC suite is a work in progress, and additional quality-of-life improvements are still being developed.

- **JavaScript for Visualizations**: To work on the underlying JavaScript for dashboard visualizations, you'll need `npx` installed. This is useful for updating or customizing the interactive features of the dashboards.

  1. **Install JavaScript Dependencies**:

     To install the necessary JavaScript packages, navigate to the directory containing the `package.json` file and run:
     ```sh
     npm install
     ```

  2. **Run Webpack**:

     To bundle the JavaScript code using Webpack, run the following command:
     ```sh
     npx webpack
     ```

---

Feel free to reach out if you need further details on PostgreSQL setup, installing `graphviz`, or configuring the HPC suite!

