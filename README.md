This project is developed to synthesize attributed social graphs that preserve 
the community structures of the input graph under differential privacy guarantee. 
For the design details of the algorithms, please refer to the paper attached. 

*******************************
Environment
*******************************
The project is implemented in JavaSE-1.8 with jgrapht-core-1.3.0. The codes are
tested on a windows machine with 8G memory. The execution does
not show any visible efficiency degradation.

*******************************
DATASETS
*******************************
Our algorithms are tested over four public datasets collected from epinions, 
facebook, petsters available on SNAP (snap.stanford.edu). For the data from 
epinions, we extract an additional dataset of 8000 nodes. Each dataset is stored
in two files: one for graph structures and the other for node attribute values.  
The corresponding file names are as follows:
   * epinions_combined.txt attribute_epinions_combined.txt
   * epinions8000_combined.txt attribute_epinions8000_combined.txt
   * facebook_combined.txt  attribtue_facebook_combined.txt  
   * petster_combined.txt attribute_petster_combined.txt

*******************************
INPUT AND OUTPUT
*******************************
The codes are used to iteratively generate synthetic attributed social graphs
based on the input dataset with 3 different methods, i.e., CAGM, DCSBM,
TriCycle. The codes guarantees differential privacy with four privacy budgets, 
i.e., 2.0, 3.0, 4.0, 5.0.   

++++++++++++++++++++++++++++++
Input parameters
++++++++++++++++++++++++++++++

parameters are required to execute the code:
    1. dataset name (string) which is selected from the following list:
        - epinions
        - epinions8000
        - facebook
        - petster
    2. starting iteration serial number (integer)
    3. ending iteration serial number (integer)
    4. method id (integer) selected from the following list:
        - 1 if CAGM will be used
        - 2 if DCSBM will be used
        - 3 if TriCycle will be used

++++++++++++++++++++++++++++++
Output files
++++++++++++++++++++++++++++++

The output of the program includes the numbered synthesised graphs and other
intermediate parameter files generated during the execution. All the output
files are stored in a folder named in the following format:

    <dataset>_[method]_combined_experiment
    
The components enclosed by angles represent the mandatory variables while those
in brackets are optional. The variable dataset is the same as the input
dataset name. The variable method is not used for CAGM and is set as TriCycle or 
DCSBM according to the method executed. This principle is also applied in the
description below.

A differentially privately synthesised graph is stored with two files:
    - the graph structure file which consists of the nodes and edges. 
    The file name TriCycle and CAGM is constructed with the following format:

        <dataset>_[method]_w_triangles_<epsilon>_final_dp_<serialno>.txt

    and the one for DCSB is in the form:

        <dataset>_DCSBM_graph_<epsilon>_final_dp_<serialno>.txt
      
    * dataset & method: the same as the naming principle of output folder name.
    * epsilon: selected from 2.0, 3.0, 4.0, 5.0.
    * serialno: a number between the starting iteration serial number and the 
      ending iteration serial number. 

    - the attribute file named in the form:
        <dataset>_[method]_attribute_<epsilon>_<serialno>.txt


*******************************
RUN THE CODES
*******************************
Please use the Jar package to run codes with the dataset files in the same
directory. The command should be in the following form: 

java -jar dpc.jar <dataset> <starting iteration serial number> <ending iteration
serial number> <method>

For example, if we want to run 10 synthetic graphs for the facebook dataset with
CAGM, the following command is used:

 java -jar dpc.jar facebook 0 9 1
 
*******************************
Evaluating measures
*******************************
All measures that evaluate the quality of synthetic graphs are implemented
in python 2.7 and stored in the file DPC_data_processing.py.
The evaluation relies on two external packages: Louvain (uploaded as
community-louvain.py) and cesna (version 4.1, Jul 25, 2018).
The Louvain package calculates the community based on only graph structures,
i.e., nodes and edges, without attributes while the cesna method takes into 
account node attributes in community detection.  



