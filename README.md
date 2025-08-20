# KINEPIK
## Overview
KINEPIK is an integrated data resource for cell signalling research. It combines multiple data sources together into one database that is programmatically accessible through an API. The data is collected both through publicly available APIs, like Omnipath, Uniprot and PubChem, but also through multiple different supplementary materials of relevant research articles. KINEPIK is publicly available at kinepik.org.

## Features
- Unified database of kinases, perturbations, phosphosites, and their interactions
- Programmatic API access
- Integration of multiple sources: OmniPath, UniProt, PubChem, and curated literature
- Easy deployment locally or on a server
- Version-controlled databases for reproducibility

## Running the app locally
To run the app locally, follow these steps:

1. Clone the repository:
`git clone https://github.com/bezzlab/KINEPIK.git`

2. Navigate to the project directory:
`cd KINEPIK`

3. Install dependencies:
`pip install -r requirements.txt`

4. Download database and version history database from [here](https://qmulprod-my.sharepoint.com/:f:/r/personal/btw695_qmul_ac_uk/Documents/__shared-QM/kinepik/database?csf=1&web=1&e=kbRvlN). Add the database to the same directory as the application. Your directory should look like this:
```
KINEPIK
├── KINEPIK_app
│   ├── api_v1
│   │   ├── __init__.py
│   │   ├── inhibitor_routes.py
│   │   ├── kinase_routes.py
│   │   ├── protein_routes.py
│   │   └── sif_routes.py
│   │
│   ├── database_scripts
│   │   ├── __init__.py
│   │   ├── API_tables.py
│   │   ├── append_phosphosites.py
│   │   ├── change_studyid.py
│   │   ├── database_structure.py
│   │   ├── discoverx_to_pert.py
│   │   ├── experimental_table.py
│   │   ├── files_tables.py
│   │   ├── gene_table.py
│   │   ├── go_terms.py
│   │   ├── kinase_type.py
│   │   ├── known_perturbations.py
│   │   ├── location_table.py
│   │   ├── perturbation_table.py
│   │   ├── phosphosite_confidence.py
│   │   ├── protein_table.py
│   │   ├── protein_table_append.py
│   │   ├── tissue_table.py
│   │   └── update_species.py
│   │
│   ├── frontend
|   |   ├── __init__.py
|   |   └── routes.py
│   │
│   ├── static
|   |   ├── css
│   │   |   └── styles.css
│   │   └── images
|   |       ├── favicon.png
|   |       ├── logo.png
│   │       └── qm-logo.png
|   |
│   ├── templates
│   │   ├── api.html
│   │   ├── index.html
│   │   └── layout.html
│   │
│   ├── __init__.py
│   └── database_con.py
│
├── testing
│   ├── json_to_excel.py
│   ├── testing_API.ipynb
│   └── testing_KINEPIK.py
│
├── application.py
├── requirements.txt
├── update_KINEPIK.py
├── KINEPIK_v0.db
└── KINEPIK_version.db
```

5. Start the application:
`python application.py`

## Example API Usage
Once the app is running locally:
- Get all kinases in the database (Database v0):  
  `http://127.0.0.1:5000/api/0/kinases/all`

- Get interactions in SIF file format at phosphosite resolution:  
  `http://127.0.0.1:5000/api/0/sif/all?resolution=phosphosite`


## Deployment
To deploy the application, follow this folder structure:
```
KINEPIK
├── KINEPIK_app
│   ├── api_v1
│   │   ├── __init__.py
│   │   ├── inhibitor_routes.py
│   │   ├── kinase_routes.py
│   │   ├── protein_routes.py
│   │   └── sif_routes.py
│   │
│   ├── database_scripts
│   │   ├── __init__.py
│   │   └── database_structure.py
│   │
│   ├── frontend
|   |   ├── __init__.py
|   |   └── routes.py
│   │
│   ├── static
|   |   ├── css
│   │   |   └── styles.css
│   │   └── images
|   |       ├── favicon.png
|   |       ├── logo.png
│   │       └── qm-logo.png
|   |
│   ├── templates
│   │   ├── api.html
│   │   ├── index.html
│   │   └── layout.html
│   │
│   ├── __init__.py
│   └── database_con.py
|
├── application.py
├── requirements.txt
├── KINEPIK_version.py
├── KINEPIK_v0.db       # preferred version of the database
└── KINEPIK_version.db
```

## Database Versioning
KINEPIK maintains a versioning system for reproducibility.  
- Default version: latest (`KINEPIK_version.db` points to the newest).  
- Access older versions by specifying them in the URL (version number):  
  `http://127.0.0.1:5000/api/0/...`

## Contributing
Contributions are welcome!  
- Fork the repo  
- Create a new branch  
- Submit a pull request with clear documentation of changes

## Citation
If you use KINEPIK in your research, please cite:
> Heinonen, V., Hübner, M., Cutillas, P., and Bessant, C. (2025). KINEPIK: An integrated data resource for cell signalling research.

## License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0).

See the [LICENSE](LICENSE) file for details.
