<div align="center">
<img src="https://github.com/hcji/ProSAP/blob/master/figure.png" width="50%">
</div>

![visitors](https://visitor-badge.glitch.me/badge?page_id=hcji/ProSAP)
![GitHub](https://img.shields.io/badge/platform-Windows%7CLinux%7CMacOS-brightgreen)
![GitHub](https://img.shields.io/github/license/hcji/ProSAP)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/hcji/ProSAP?include_prereleases)
![GitHub top language](https://img.shields.io/github/languages/top/hcji/ProSAP)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5763315.svg)](https://doi.org/10.5281/zenodo.5763315)

ProSAP (Protein Stability Analysis Pod) is standalone and user-friendly software with graphical user interface (GUI). 
ProSAP provides an integrated analysis workflow for thermal shift assay, which includes five 
modules: data preprocessing, data visualization, TPP analysis, NPARC analysis and iTSA analysis. 
With the assistance of the user-friendly interface, researchers can easily compare several statistical 
strategies, analyze the results and draw the conclusion from the proteomics quantitative table obtained 
by Proteome Discoverer or MaxQuant. Users would also benefit from a comprehensive overview of the 
performance of different algorithms, and apply appropriate algorithms to their dataset easily.


## Installation

**Stable version:**

Windows: [ProSAP-0.99.5b-Windows.zip](https://zenodo.org/record/5763315/files/ProSAP-0.99.5b-Windows.zip)    
Linux: [ProSAP-0.99.5b-Linux.tar.xz](https://zenodo.org/record/5763315/files/ProSAP-0.99.5b-Linux.tar.xz)    
MacOS: [ProSAP-0.99.5b-Mac.zip](https://zenodo.org/record/5763315/files/ProSAP-0.99.5b-Mac.zip)

#### For Windows: 
The installer version should be preferred but might require administrator permissions. 
Please follow the video of the installation at [video folder](https://github.com/hcji/ProSAP/tree/master/video).

**Note**: 
1. Since we do not pay Microsoft for certification, you might have to confirm that you want to trust 
"software from an unknown source". 
2. When installing, please do **not** change the default path, because sometimes we find it does not run correctly with customized path.  
3. We observed unknown errors happened if you have installed multiple versions of R. Uninstall all version of R and reinstall the latest version will fix the error.  
4. When running for the first time, it may need a long time to install the R package. There will be some requirements showing in the command prompt window, you should just follow the instructions.  

#### For Linux:  
Please follow the following installation steps and refer the video of the installation at [video folder](https://github.com/hcji/ProSAP/tree/master/video):   

1. Install [R](https://cran.r-project.org/)

        sudo apt-get install r-base
  
2. Install extra dependency of QT:

        sudo apt-get install libxcb-xinerama0

3. Download the tgz file and unzip to any folder and execute ProSAP.  
4. Execute the binary file:

        ./ProSAP
        
5. When running for the first time, it may need a long time to install the R package. There will be some requirements showing in the command prompt window, you should just follow the instructions.  
        
        
#### For MacOS:  
Please follow the following installation steps and refer the video of the installation at [video folder](https://github.com/hcji/ProSAP/tree/master/video): 

1. Allow software from unknown developers:  

        sudo spctl --master-disable
        
2. Install [R](https://cran.r-project.org/)
3. Download the zip file and unzip to any folder and execute ProSAP.
4. Execute the binary file:

        ./ProSAP
        
5. When running for the first time, it may need a long time to install the R package. There will be some requirements showing in the command prompt window, you should just follow the instructions.  


#### For MacOS with M1 chip:
1. If you are using MacBook Pro with M1 chip, you should install [Rosetta](https://machow2.com/rosetta-mac/).    

        softwareupdate --install-rosetta --agree-to-license

2. Install [R](https://www.r-project.org/) of **High Sierra version** (Inter 64 based). This is because rpy2 has not support ARM64 based R yet. 
We will update if rpy2 support it in the future.  
3. Follow the same step 1 ~ step 5 of MacOS.  


**Note:** ProSAP has been test on Windows 7, Windows 10, Windows 11, Ubuntu 20.04 and MacOS 11.6.1. However, it does not work on Windows XP.

## Development version

#### For Windows or MacOS::  
1. Install [Anaconda](https://www.anaconda.com/)  or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)   
2. Install [R](https://cran.r-project.org/)  
3. Install [Git](https://git-scm.com/downloads)  
4. Open commond line, create environment and enter with the following commands:  

        conda create -n ProSAP python=3.8
        conda activate ProSAP

**Note:** MacBook Pro with **M1** chip may need python=3.9 and make sure R version is suitable.  
5. Clone the repository and enter:  

        git clone https://github.com/hcji/ProSAP.git
        cd ProSAP

6. Install dependency with the following commands:  
        
        pip install requirements.txt
        
7. Run ProSAP.py:  

        python ProSAP.py
        
8. Or, you may want to compile binary and run exe:  

        pyinstaller ProSAP.py -i ./img/ProSAP.ico --hidden-import=“sklearn.utils._weight_vector” 
        cd dist/ProSAP
        ProSAP.exe
        
#### For Linux: 
1. Most Linux distributions have included git and conda, but you may need install extra dependency of QT:  

        sudo apt-get install libxcb-xinerama0
        
2. Install [R](https://cran.r-project.org/)

        sudo apt-get install r-base
        
3. Refer the step 3 ~ step 7 of Windows or MacOS.  


## Document

The detailed usage is included in the [ProSAP website](https://hcji.github.io/ProSAP_Pages/).    
The videos for using the software are available at the [video folder](https://github.com/hcji/ProSAP/tree/master/video).    

## Contact

Ji Hongchao   
E-mail: ji.hongchao@foxmail.com    

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0002-7364-0741" href="https://orcid.org/0000-0002-7364-0741" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0002-7364-0741</a></div>
