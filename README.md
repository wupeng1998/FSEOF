



# Installation & Usage

1. Make sure you have python installed on your computer
2. Download the FSEOF GitLab repository as a zip file
3. Unzip the downloaded zip file
4. Move the "FSEOF" folder to the desired destination on your computer
5. Paste the .xml file into the folder "FSEOF"
6. Navigate to the folder in the terminal of your computer. Users without knowledge of navigation in the terminal have the following options:
7. In the now opened terminal window type in the following line and press enter. **This is only necessary for the first time you use the tool!**
    
    `pip install -r requirements.txt`
  
8. Type in `python run_FSEOF.py NameOfYourSBMLFile BiomassID reactionID` and press enter. The results are stored in an Excel file
> **Example:** `python run_FSEOF.py yeast_gem.xml r_4041 r_426`

> **OR:** python FSEOF_self.py yeast_gem.xml biomass r_4041
