# ProSAP

![GitHub](https://img.shields.io/badge/platform-Windows%7CLinux%7CMacOS-brightgreen)
![GitHub](https://img.shields.io/github/license/hcji/ProSAP)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/hcji/ProSAP?include_prereleases)
![GitHub top language](https://img.shields.io/github/languages/top/hcji/ProSAP)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5525802.svg)](https://doi.org/10.5281/zenodo.5525802)


ProSAP (Protein Stability Analysis Pod) is standalone and user-friendly software with graphical user interface (GUI). 
ProSAP provides an integrated analysis workflow for thermal shift assay, which includes five 
modules: data preprocessing, data visualization, TPP analysis, NPARC analysis and iTSA analysis. 
With the assistance of the user-friendly interface, researchers can easily compare several statistical 
strategies, analyze the results and draw the conclusion from the proteomics quantitative table obtained 
by Proteome Discoverer or MaxQuant. Users would also benefit from a comprehensive overview of the 
performance of different algorithms, and apply appropriate algorithms to their dataset easily.

<div align="center">
<img src="https://github.com/hcji/ProSAP/blob/master/figure.png" width="50%">
</div>


## Installation

Download links:

Windows: [ProSAP-0.99.3-win64.zip](https://zenodo.org/record/5525802/files/ProSAP-0.99.3-win64.zip?download=1)    
Linux: [ProSAP-0.99.3-Linux.tar.gz](https://zenodo.org/record/5525802/files/ProSAP-0.99.3-Linux.tar.gz?download=1)    
MacOS: [ProSAP-0.99.3-Mac.zip](https://zenodo.org/record/5525802/files/ProSAP-0.99.3-Mac.zip?download=1)

For Windows, the installer version should be preferred but might require administrator permissions. 
Since we do not pay Microsoft for certification, you might have to confirm that you want to trust 
"software from an unknown source". When installing, please do **not** change the default path.   

**Note**: We observed unknown errors happened if you have installed multiple versions of R. Uninstall all version of R
and reinstall the latest version will fix the error. 


For Linux, you should install [R](https://cran.r-project.org/) first, Then, download the tgz file, 
unzip to any folder and execute ProSAP. You may also need install extra dependency of QT:

        sudo apt-get install libxcb-xinerama0

For MacOS, you should install [R](https://cran.r-project.org/) first, Then, download the zip file, 
unzip to any folder and execute ProSAP. Since we do not pay Apple for certification, 
you might have to allow software from unknown developers:

        sudo spctl --master-disable


ProSAP has been test on Windows 7, Windows 10, Windows 11, Ubuntu 20.04 and MacOS 10.15.7. However, it does not work on Windows XP.

## Compile with source

Advanced users can compile the source codes. Required dependencies:

* [Anaconda for python (python version >= 3.6)](https://www.anaconda.com/)    
* [R (>= 4.0)](https://www.r-project.org/)    
* [PyQt5](https://pypi.org/project/PyQt5/)    
* [RPy2](https://pypi.org/project/rpy2/)
* [pyinstaller](https://www.pyinstaller.org/)    
* [mkl](https://pypi.org/project/mkl/) **(for MacOS only)**

Then, clone the repository and enter:
        
        git clone https://github.com/hcji/ProSAP.git

Next, compile the source   
        
        pyinstaller ProSAP.py -i ./img/ProSAP.ico --hidden-import=“sklearn.utils._weight_vector” 

Just run ./dist/ProSAP


## Document

The detailed usage is included in the [ProSAP website](https://hcji.shinyapps.io/prosap_page/).    
The videos for using the software are available at the [video folder](https://github.com/hcji/ProSAP/tree/master/video).    

## Contact

Ji Hongchao ![Twitter URL](https://img.shields.io/twitter/url?url=https%3A%2F%2Ftwitter.com%2FHongchaoJ%2Fstatus%2F1440875003478564866)    
E-mail: ji.hongchao@foxmail.com   

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0002-7364-0741" href="https://orcid.org/0000-0002-7364-0741" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0002-7364-0741</a></div>