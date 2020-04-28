# Development of a data science machine learning platform for food quality
**Background** : Approximately 88 million tonnes of food are wasted annually across the European Union due to spoilage, incurring serious financial and environmental implications. The safety of poultry-containing food products can be compromised when the meat is prepared and packaged. The presence of pseudomonads in aerobic conditions can produce unpleasant slime and off-odours compromising the quality of the meat. European Commission guidelines dictate that the quality of fresh meat be determined by total viable counts (TVC) of bacteria or counts of Enterobacteriaceae. Obtaining this type of data for meat samples is time-consuming and labor intensive. Multispectral imaging (MSI) and Fourier Transform Infrared (FTIR) spectroscopy technologies are being developed to provide real-time, non-invasive, on-line screening for microbiological safety and spoilage assessment of meat samples. 

**Aims and objectives** : To create a scalable platform based on Representational State Transfer (REST) and application of APIs (Application Programming Interface) for the rapid quality assessment of chicken-thigh and chicken-burger samples. This will be achieved through the implementation of machine learning regression algorithms based on MSI and FTIR data. 

**Approach** : The approach includes the employment of multiple machine learning algorithms to achieve the optimum prediction accuracies. In particular, five algorithms have been implemented: Linear Regression, K-Nearest Neighbor, Random Forest, Support Vector Machines with Polynomial Kernel and Radial Basis Function Kernel. Analyses and modelling were conducted via R language and environment for statistical computing and graphics, including the API (Plumber) which acts as a “pipeline” facilitating data requests and returning a response as a predicted value of the bacterial count present in a given sample. 

## Table of Contents
* [Fresh-API](#Fresh-API)
* [Installation](#Installation)
  * [Dependencies](#Dependencies)
  * [Install Fresh-API from source](#Install Fresh-API from source)
  * [Step 1 — Installing Docker](#Step 1 — Installing Docker)
  * [Step 2 — Building the Docker Image](#Step 2 — Building the Docker Image)
  * [Step 3 — Running the Docker container](#Step 3 — Running the Docker container)
  * [Step 4 — Sending request to REST-API](#Step 4 — Sending request to REST-API)
* [JSON file format specification](#JSON file format specification)
* [Additional commands for Docker](#Additional commands for Docker)



## Fresh-API 
The purpose of the application is to provide a platform to predict the total bacterial counts present on a given meat sample through the implementation of REST API. Fresh-API is build upon the technology of RESTful web based API consisting of an **Endpoint URL** and **Body Data** which can be transmitted via **HTTP request methods** by sending a single JSON-encoded data string. Machine Learning regression models are developed onto the platform and validated against new data coming to the system by sending the prediction results through the REST backend.

# Installation
## **Dependencies**
 * UNIX Operating Systems  
 * [Docker](https://www.docker.com/why-docker)
 * R (tested on version 3.6.3)
 * The following R libraries:
 ```
 
 ```
 
## **Install Fresh-API from source**  
```
```

## **Step 1 — Installing Docker**

To download the latest version, install Docker from the official Docker repository. This section guides you how to do that.

First, in order to ensure the downloads are valid, add the GPG key for the official Docker repository to your system:
```
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```
Add the Docker repository to APT sources:
```
$ sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```
Update the package database with the Docker packages from the newly added repository:
```
$ sudo apt-get update
```
To ensures that you are running the installation from the official Docker repository, run the command:
```
$ apt-cache policy docker-ce
```
This should give an output similar to the following:
```
docker-ce:
  Installed: 5:19.03.8~3-0~ubuntu-xenial
  Candidate: 5:19.03.8~3-0~ubuntu-xenial
  Version table:
 *** 5:19.03.8~3-0~ubuntu-xenial 500
        500 https://download.docker.com/linux/ubuntu xenial/stable amd64 Packages
        100 /var/lib/dpkg/status
     5:19.03.7~3-0~ubuntu-xenial 500
        500 https://download.docker.com/linux/ubuntu xenial/stable amd64 Packages
 ```
From the output, you will notice that the docker-ce is not yet installed. However, the output will show the target operating system and the version number of the Docker. Please note that version numbers may differ depending on the time of installation.

Use the following command to install Docker:
```
sudo apt install docker-ce
```
This will install Docker, start the daemon and enable it to automatically start on boot. To confirm that the Docker is active and working, run:
```
sudo systemctl status docker
```
The command will provide the following output:
```
● docker.service - Docker Application Container Engine
   Loaded: loaded (/lib/systemd/system/docker.service; enabled; vendor preset: enabled)
   Active: active (running) since Mon 2020-03-23 12:34:27 GMT; 4 weeks 1 days ago
     Docs: https://docs.docker.com
 Main PID: 11575 (dockerd)
    Tasks: 31
   Memory: 1.9G
      CPU: 33min 41.448s
   CGroup: /system.slice/docker.service
           └─11575 /usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock
```           
## **Step 2 — Building the Docker Image**
`docker` consists of chain of options and commands followed by arguments. The syntax takes this form:
```
$ docker [option] [command] [arguments]
```
To view all available subcommands, type:
```
$ docker
```
The `docker build` command builds Docker images from a Dockerfile and a “context”. A build’s context is the set of files located in the specified `PATH`:
```
$ docker build -t freshapi .
```
NOTE: `-t` or `--tag` is a `docker build` option used for label a docker image in the `name:tag` format. In this case, name will be `freshapi`. The `.` operator specifies the current directory containing the 'Dockerfile'. Optionally, the version of an image can be tagged in which the user would like to run the container. For example, `docker build -t freshapi:2.0 .` 

## **Step 3 — Running the Docker container**
The `docker run` command creates a container from a given image and starts the container using a given command or entrypoint. To run the resulting docker image, the command is as follows:
```
$ docker run --name [container name of choice] -d freshapi
```
NOTE: `-d` or `--detach` is a `docker run` option for running the container in background. `--name` is used to assign a name to the container. 

In order to execute commands on running containers, `docker exec` is used to specify the container name as well as the command to be executed on this container:
```
$ docker exec -it [container name of choice] bash
```
NOTE: The combination of `-i` and `-t` allows the container to run in an interative mode and provides access to the terminal to execute further commands in the application. The `bash` command creates a new Bash session for the container.

## **Step 4 — Sending request to REST-API**
The design and features of Fresh-API is to accept request inputs as `.json` data and query parameters defined by the developers. `curl` is a command line tool for transferring and receiving HTTP requests and responses. The syntax for the four endpoints are as follows:
1. Bacterial count prediction of MSI data:
```
curl --data [MSI(.json)] "localhost:80/predict_MSI?bacteria=&model="
```
2. Bacterial count prediction of FTIR data:
```
curl --data [FTIR(.json)] "localhost:80/predict_FTIR?product=&bacteria=&model="
```
3. Bacterial count prediction using the most accurate Machine Learning method:
```
curl --data [MSI or FTIR(.json)] "localhost:80/predict?platform=&product="
```
4. Report of Machine Learning models for MSI and FTIR data:
```
curl "localhost:80/report?platform=&product="
```
The following table provides a list of endpoints and their corresponding query parameters:
|Endpoints|Query parameters| 
|---------|----------------|
|`#* @post /predict_MSI`|`bacteria= TVC or Ps` ; `model= rf, knn, svmLinear, svmRadial, svmPoly, lm`|
|`#* @post /predict_FTIR`|`product= CB, CTF` ; `bacteria= TVC, Ps, Bth, LAB` ; `model= rf, knn, svmLinear, svmRadial, svmPoly, lm`|
|`#* @post /predict`|`platform= MSI, FTIR` ; `product= CB, CTF`|
|`#* @get /report`|`platform= MSI, FTIR` ; `product= CB, CTF`|

NOTE: `--data` or `-d` denotes the `curl` command used for passing data to the request body. `platform`, `product`, `bacteria` and `model` parameters were passed to the endpoint using “query strings”. The (`?`) appended to the URL indicates the start of the query string. In the query string, each parameter is concatenated with other parameters through the ampersand (`&`) symbol.
# JSON file format specification
* In the case of predicting bacterial counts, data derived from an analytical platform such as MSI should contain the 18 mean values as features at the beginning of a JSON file.
* For FTIR derived samples the file should contain wavelengths in the range of 1001-4000 nm. 
# Additional commands for Docker
| USAGE | Commands |
|-------|----------|
|Display logs of a container|`$ docker logs [container name]`|
|List all existing containers (running and not running)|`$ docker ps -a`|
|List your images|`$ docker image ls`|
|Stop a specific container|`$ docker stop [container name]`|
|Stop all running containers|`$ docker stop $(docker ps -a -q)`|
|Delete a specific container (only if stopped)|`$ docker rm [container name]`|
|Delete all containers (only if stopped)|`$ docker rm $(docker ps -a -q)`|
|Delete all unused containers, unused images and networks|`$ docker system prune -a --volumes`|


