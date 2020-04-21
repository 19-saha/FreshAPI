# Development of a data science machine learning platform for food quality
**Background** : Approximately 88 million tonnes of food are wasted annually across the European Union due to spoilage, incurring serious financial and environmental implications. The safety of poultry-containing food products can be compromised when the meat is prepared and packaged. The presence of pseudomonads in aerobic conditions can produce unpleasant slime and off-odours compromising the quality of the meat. European Commission guidelines dictate that the quality of fresh meat be determined by total viable counts (TVC) of bacteria or counts of Enterobacteriaceae. Obtaining this type of data for meat samples is time-consuming and labor intensive. Multispectral imaging (MSI) and Fourier Transform Infrared (FTIR) spectroscopy technologies are being developed to provide real-time, non-invasive, on-line screening for microbiological safety and spoilage assessment of meat samples. 

**Aims and objectives** : To create a scalable platform based on Representational State Transfer (REST) and application of APIs (Application Programming Interface) for the rapid quality assessment of chicken-thigh and -burger samples. This will be achieved through the implementation of machine learning regression algorithms based on MSI and FTIR data. 

**Approach** : The approach includes the employment of multiple machine learning algorithms to achieve the optimum prediction accuracies. In particular, five algorithms have been implemented: Linear Regression, K-Nearest Neighbor, Random Forest, Support Vector Machines with Polynomial Kernel and Radial Basis Function Kernel. Analyses and modelling were conducted via R language and environment for statistical computing and graphics, including the API (Plumber) which acts as a “pipeline” facilitating data requests and returning a response as a predicted value of the bacterial count present in a given sample. 

## Table of Contents
* [Food Spoilage-API ](#Food-Spoilage-API)




## Food Spoilage-API 
The purpose of the application is to provide a platform to predict the total bacterial counts present on a given meat sample through the implementation of REST API. Food Spoilage-API is build upon the technology of RESTful web based API consisting of an **Endpoint URL** and **Body Data** which can be transmitted via **HTTP request methods** by sending a single JSON-encoded data string. Machine Learning regression models are developed onto the platform and validated against new data coming to the system by sending the prediction results through the REST backend.

