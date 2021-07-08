## TSAnalyst

TSAnalyst is standalone and user-friendly software with graphical user interface (GUI). TSAnalyst 
TSAnalyst provides an integrated analysis workflow for thermal shift assay, which includes five 
modules: data preprocessing, data visualization, TPP analysis, NPARC analysis and iTSA analysis. 
With the assistance of the user-friendly interface, researchers can easily compare several statistical 
strategies, analyze the results and draw the conclusion from the proteomics quantitative table obtained 
by Proteome Discoverer or MaxQuant. Users would also benefit from a comprehensive overview of the 
performance of different algorithms, and apply appropriate algorithms to their dataset easily.

<div align="center">
<img src="https://github.com/hcji/TSAnalyst/blob/master/figure.png">
</div>


### Installation

For Windows, the installer version should be preferred but might require administrator permissions. 
Since we do not pay Microsoft for certification, you might have to confirm that you want to trust 
"software from an unknown source". For Linux, you should install [R](https://cran.r-project.org/) first, 
Then, download the gzip file, unzip to any folder and execute TSAnalyst.

Windows: [TSAnalyst-0.99.2-win64.exe](https://zenodo.org/record/5077939#.YOZHmDPitPY)    
Linux: [TSAnalyst-0.99.2-Linux.tar.gz](https://zenodo.org/record/5077939#.YOZHmDPitPY)

### Compile with source

Advanced users and Mac users can compile the source codes. Required dependencies:

* [Anaconda for python (python version >= 3.6)](https://www.anaconda.com/)    
* [R (>= 4.0)](https://www.r-project.org/)    
* [PyQt5](https://pypi.org/project/PyQt5/)    
* [RPy2](https://pypi.org/project/rpy2/)    

Then, clone the repository and enter:
        
        git clone https://github.com/hcji/TSAnalyst.git

Next, compile the source   
        
        pyinstaller TSAnalyst.py -w -i ./img/TPCA.ico --hidden-import=“sklearn.utils._weight_vector” 

Just run */dist/TSAnalyst.exe*


### Document

The detailed usage is included in the [help document](https://github.com/hcji/TSAnalyst/blob/master/TSAnalyst%20Desktop%20Tutorial.docx)


### Contact

E-mail: ji.hongchao@foxmail.com
