# tcga-data-nf
## End-to-end reproducible processing of cancer regulatory networks from TCGA data

![](https://github.com/QuackenbushLab/tcga-data-nf/workflows/build/badge.svg)

version = '0.0.17'

Workflow to download and prepare TCGA data.

The workflow divides the process of downloading the data in two steps:
1. Downloading the raw data from GDC and saving the rds/tables needed later
2. Preparing the data. This step includes filtering the data, normalizing it... 
3. Analysis of gene regulatory networks

The idea is that data should be downloaded once, and then prepared for the task at hand.

>**Where can I find more details about this workflow?**  
If you want more details about the scope and use of this workflow, for instance, you want to decide if it's useful for your research, we
recommend you check out the paper: ["Reproducible processing of TCGA regulatory
networks"](https://www.biorxiv.org/content/early/2024/11/07/2024.11.05.622163).

>**Are there examples of how to configure the workflow or sample datasets?**  
Of course! In the [QuackenbushLab/tcga-data-supplement](https://github.com/QuackenbushLab/tcga-data-supplement)
repository you will find the companion data and configuration files for the paper. You can read about a full analysis we
did on colon cancer and find all links/instructions for the precomputed GRNs of common cancer types.

:warning: GitHub is still not allowing the QuackenbushLab to run actions. The docker image is built and pushed 
manually.


## Getting started 

1. First you'll need to [install](https://www.nextflow.io/docs/latest/install.html) nextflow on your machine. Follow the
   `hello world` example to check if Nextflow is up and running.
2. Pull the workflow `nextflow pull QuackenbushLab/tcga-data-nf`
3. Install and pull the docker/singularity container or conda to run the whole pipeline  
  3a. **DOCKER**: 
    - [Install](https://docs.docker.com/engine/install/) docker on your machine
    - Pull the container for this workflow `docker pull violafanfani/tcga-data-nf`
    - use the `-profile docker` option when running the workflow  
  3b. **CONDA**:    
    - Install [conda](https://docs.anaconda.com/miniconda/).
    - Use the `-profile conda` option when running the workflow
4. Run some test workflows  
  4a. test the download: `nextflow run QuackenbushLab/tcga-data-nf -profile <docker/conda>,testDownload `  
  4a. test the prepare: `nextflow run QuackenbushLab/tcga-data-nf -profile <docker/conda>,testPrepare`  
  4a. test the analyze: `nextflow run QuackenbushLab/tcga-data-nf -profile <docker/conda>,testAnalyze`  
  4a. test the full workflow: `nextflow run QuackenbushLab/tcga-data-nf -profile <docker/conda>,test `  

If you can run all these steps, you can procede defining your own configuration files and run your own analysis. 

:warning: The full workflow can be slow (>45minutes).  

Check the docs for [AWS](docs.md#AWS) for the steps on how to run the workflow on a simple EC2 instance.
These steps could also help as a quickstart to check that you have everything up and running.

## Running the workflow

### Docker

The docker container is hosted on docker.io. 

```
docker pull violafanfani/tcga-data-nf:0.0.17
```

More details on the container can be found in the [docs](docs.md#Docker)

### Conda

Alternatively, one can run the workflow with conda environments. 
In order to create and use conda one can pass it as a profile `-profile conda`
as:

```
nextflow run QuackenbushLab/tcga-data-nf -profile conda,test ...
```

For the moment we are using one single environment to be used with all the r
scripts. This allows the pipeline to generate the environment only once (which can 
be time consuming) and then to reuse it.

The three environment are inside the `containers/conda_envs` folder: 
1. merge_tables
2. r_all
3. analysis

However, to improve portability, we use the process selectors labels to specify the different environments, allowing the user to specify their own environments too. 

More details in the [docs](docs.md#conda).


### CPU and memory requirement

CPU and memory requirements depend on the size of the datasets, but we report a whole execution report for the full
pipeline. 

You can find that at [full](execution_report_test_full.html) execution report that was obtained by running the v0.0.17
version of tcga-data-nf with standard configuration files and no parallelization of processes on AWS EC2 instance c5.4xlarge with 32Gb of memory and 16 vCPUs. 

### Install or update the workflow

Before running the workflow we recommend pulling the last version with the following command

```bash
nextflow pull QuackenbushLab/tcga-data-nf
```

### Run the analysis

One can run the workflow by simply using the `nextflow run` command and using a custom configuration file.

```bash
nextflow run QuackenbushLab/tcga-data-nf -c my-config.conf
```

Below we give more details on the configuration steps

## Configuration

First there are three main parameters that need to be passed by the user:
- `resultsDir = "results"`: general folder under which you want to find the results. This can directly reference an AWS S3 bucket
- `batchName = "my-batch"`: name of the run, this is gonna create a subfolder where the results are stored.
- `pipeline = 'download'`: name of the pipeline, one of download,prepare,analyze,full

This way all data generated by the pipeline will be found inside the `resultsDir/batchName/` folder.
If nothing is passed, all results will be in the `results/my-batch` folder.

Secondly, you'll need a 


For a full list of the configuration parameters check [here](docs.md##configurations).

## Results

Below is the structure you can expect in the output folder when you run each pipeline. 

For each case we report the output of the testing profile:
- download pipeline: `-profile testDownload`
- prepare pipeline: `-profile testPrepare`
- analyze pipeline: `-profile testAnalyze`
- full pipeline: `-profile test`

Detailed output folder structure can be found at the [docs](docs.md##result-folders)

## Development

In case you wanted to make modifications to the workflow and/or run it locally

0. Fork the repo into your own github
  
1. Clone the forked nextflow repo 
   ```bash
   git clone git@github.com:myaccount/tcga-data-nf.git
   ```
2. Build docker locally 
    ```bash
    docker build . -f containers/Dockerfile --build-arg CONDA_FILE=containers/env.base.python.yml --no-cache -t my-tcga-data-nf:latest
    ```
3. In alternative, use the conda profile
4. Run your workflow
    ```bash
    nextflow run .  -profile testDownload --resultsDir myresults/ --pipeline download -with-docker my-tcga-data-nf:latest
    ```


### Add a tool

[extend_workflow](extend_workflow.md)

## Authors

Maintainer of the workflow:
- Viola Fanfani (vfanfani@hsph.harvard.edu)

Maintainer of NetworkDataCompanion:
- Kate Hoff Shutta

Other contributors:
- Panagiotis Mandros
- Jonas Fischer
- Soel Micheletti
- Enakshi Saha
- Chen Chen

## Cite

The companion preprint is now on bioRxiv:

Reproducible processing of TCGA regulatory networks
Viola Fanfani, Katherine H. Shutta, Panagiotis Mandros, Jonas Fischer, Enakshi Saha, Soel Micheletti, Chen Chen, Marouen Ben Guebila, Camila M. Lopes-Ramos, John Quackenbush
bioRxiv 2024.11.05.622163; doi: https://doi.org/10.1101/2024.11.05.622163

```latex
@article {Fanfani2024.11.05.622163,
	author = {Fanfani, Viola and Shutta, Katherine H. and Mandros, Panagiotis and Fischer, Jonas and Saha, Enakshi and Micheletti, Soel and Chen, Chen and Ben Guebila, Marouen and Lopes-Ramos, Camila M. and Quackenbush, John},
	title = {Reproducible processing of TCGA regulatory networks},
	elocation-id = {2024.11.05.622163},
	year = {2024},
	doi = {10.1101/2024.11.05.622163},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/11/07/2024.11.05.622163},
	eprint = {https://www.biorxiv.org/content/early/2024/11/07/2024.11.05.622163.full.pdf},
	journal = {bioRxiv}
}
```

## Changelog

### v0.0.17

:warning: DRAGON :dragon:: Be careful, now DRAGON allows to run it with any input. By default, in the full pipeline we run it for
both methylation-expression and CNV-expression. Here are the things you'll find different: 
- Input: the configuration metadata will now need 5 fields `uuid,type1,methylation,type2,expression`. This way the used
  can now specify what is the type of the two files, for instance methylation-expression or cnv-expression. 
  Also, if type2 is not "expression" dragon will assume you already aligned the data and will run it with any data type
  as long as they are compatible tables.
- Output: the files will now have structure `uuid_type1_type2_dragon`


:dromedary_camel: ALPACA: The ALPACA method is now available and can be run just by specifying the right parameters in the
configuration file

:sparkles: Other GRN: We have now added GENIE3 and WGCNA to the methods, and included a #howto guide on 
[how to add new methods](extend_workflow.md).