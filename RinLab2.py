#Description: This code is made to extract data from the .OUT files from SIMULATE 3D and save the results as a .MAT file
# the data is stored as a 3D numpy arrays with dimensions (no_radial_nodes,, no_radial_nodes, no_axial_nodes)
# the data is saved as MAT files using scipy.io.savemat
from distutils.log import error
import os
import re
from tkinter import filedialog 
from tkinter import *
import glob
import numpy as np
import scipy.io as sio


# Global variables

import re

def initialise_global_variables(filename:str):
    """
    Initializes global variables based on the contents of a file.

    Parameters:
    filename (str): The path to the file containing the data.

    Returns:
    None
    """
    global no_radial_nodes
    global lines
    global no_axial_nodes
    global no_cases
    global axial_nodes
    global radial_nodes
    global length
    global cases
    global n_timesteps
    global list_length
    global time_array
    global read_test_FLX1
    global Fast_flux_changed

    with open(filename, "rt") as myfile:
        # read in the string and find the number of radial nodes    
        lines = myfile.readlines() # read in the file as a list of lines
    try:
        no_radial_nodes = int(re.findall(r'QTR.MAP Dimensions . . . . . . .ID      (\d\d)+', ''.join(lines))[0])
    except IndexError:
        no_radial_nodes = 34
        print(f"Number of radial nodes not found assuming {no_radial_nodes} radial nodes")
        pass       
    try:
        no_axial_nodes = int(re.findall(r'(?<=KD\s\s)(\d\d)', ''.join(lines))[0])
    except IndexError:
        no_axial_nodes = 26
        print(f"Number of axial nodes not found assuming {no_axial_nodes} axial nodes")
        pass
    try:
        no_cases = len(set(re.findall(r'Case  [1-9]+', ''.join(lines)))) - 1
    except IndexError:
        error("No case numbers found")

    axial_nodes = range(1, no_axial_nodes + 1) # range of axial nodes with matlab indexing
    radial_nodes = range(1, no_radial_nodes + 1) # range of radial nodes with matlab indexing
    cases = range(1, no_cases + 1) # range of cases with matlab indexing
    length = len(lines) # number of lines in the file
    t0 = 0 
    tf = 75 #hours
    dt = 0.5 #hours
    n_timesteps = int((tf-t0)/dt)
    list_length = int(no_radial_nodes*no_radial_nodes*no_axial_nodes)
    time_array = np.arange(t0,tf+dt,dt)
    read_test_FLX1 = TRUE
    Fast_flux_changed = FALSE
            
# Define globals methods

# Get filename functions

import glob

def is_only_one_out_file_in_current_directory():
    """
    Check if there is only one '.OUT' file in the current directory.

    Returns:
        bool: True if there is only one '.OUT' file, False otherwise.
    """
    return len(glob.glob('*.OUT')) == 1
def get_file_path():
    """
    Opens a file dialog to select a file and returns the selected file path.

    Returns:
        str: The selected file path.
    """
    root = Tk()
    root.withdraw()
    filename = filedialog.askopenfilename(initialdir = os.getcwd(),title = "Select file",filetypes = (("OUT files","*.out"),("all files","*.*")))
    root.destroy()
    return filename
def get_out_file_path():
    if is_only_one_out_file_in_current_directory():
        return glob.glob('*.OUT')[0]
    else:
        return get_file_path()
    
# Get power level
import re

def get_power_level(line: str) -> float:
    """
    Extracts the power level from a given line.

    Args:
        line (str): The line containing the power level information.

    Returns:
        float: The extracted power level.

    Raises:
        IndexError: If no power level is found in the line.
    """
    return float(re.findall(r'CTP\s+(\d+.\d+)', ''.join(line))[0])
def get_relative_power_level(line:str):
    """
    Extracts the relative power level from a given line.

    Args:
        line (str): The line containing the power level information.

    Returns:
        float: The relative power level.

    Raises:
        IndexError: If the power level cannot be found in the line.
    """
    return float(re.findall(r'PERCTP\s+(\d+.\d+)',''.join(line))[0])
# Get average exposure
def get_average_exposure(line:str):
    """
    Calculates the average exposure from a given line.

    Parameters:
    line (str): A string containing the line of data.

    Returns:
    float: The average exposure value extracted from the line.

    Raises:
    IndexError: If no exposure value is found in the line.

    """
    return float(re.findall(r'EBAR\s+(\d+.\d+)',''.join(line))[0])
# verify case line
def case_line_identifier(line:str):
    """
    Checks if a given line contains a case number.

    Parameters:
    line (str): A string containing the line of data.

    Returns:
    bool: True if the line contains a case number, False otherwise.

    """
    identifier = re.search(r'Case\s+\d+\s+Step\s+\d+\s+Ringhals 4\s+cycle\s+\d+',line)
    if identifier:
        return True
    else:
        return False
def variable_line_identifier(line:str,variable_match):
        """
        Checks if a given line constains a variable name.

        Parameters:
        line (str): A string containing the line of data.

        Returns:
        bool: True if the line constains a variable name, False otherwise.

        """
        if type(variable_match) == re.Pattern:
            identifier = variable_match
        elif type(variable_match)==str:
            identifier = re.escape(variable_match)
        else:
            error("Passed variable_match is not a string or a regular expression")
        if re.search(identifier,line):
            return True
        else:
            return False
        
def Case_search_activate(): #Acitvate the search for the case number and deactivate other searches
    """
    Activates the search for the case number and deactivates other searches.
    
    Returns:
    None
    """
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary
    search_case = True
    search_variable = False
    search_indices = False
    search_data = False
    search_summary = False

def Variable_search_activate(): #Acitvate the search for the variable and deactivate other searches
    """
    Activates the search for the variable and deactivates other searches.
    
    Returns:
    None
    
    """
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary

    search_case = False
    search_variable = True
    search_indices = False
    search_data = False
    search_summary = False
def Indices_search_activate(): #Acitvate the search for the indices and deactivate other searches
    """
    Activates the search for the indices and deactivates other searches.
    
    Returns:
    None
    
    """
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary
    search_case = False
    search_variable = False
    search_indices = True
    search_data = False
    search_summary = False
def Data_search_activate(): #Acitvate the search for the data and deactivate other searches
    """
    Activates the search for the data and deactivates other searches.
    
    Returns:
    None
    
    """
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary
    search_case = False
    search_variable = False
    search_indices = False
    search_data = True
    search_summary = False
def Step_search_activate(): #Acitvate the search for the data and deactivate other searches
    """
    Activates the search for the data and deactivates other searches.
    
    Returns:
    None
    
    """
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary
    search_case = False
    search_variable = False
    search_indices = False
    search_data = False
    search_summary = False
def Summary_search_activate():
    global search_case
    global search_variable
    global search_indices
    global search_data
    global search_summary
    search_case = False
    search_variable = False
    search_indices = False
    search_data = False
    search_summary = True
# Get case number
def get_case_number(line:str): # Find the case number in the beginning of the line ignore later case numbers
    """
    Extracts the case number from a given line.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    int: The extracted case number.
    
    Raises:
    IndexError: If no case number is found in the line.
    
    """
    return int(re.findall(r'Case\s+(\d+)',line)[0])
def get_step_number(line:str): # Find the step number 
    """
    Extracts the step number from a given line.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    int: The extracted step number.
    
    Raises:
    IndexError: If no step number is found in the line.
    
    """
    return int(re.findall(r'Step\s+(\d+)',line)[0])
def get_time(line:str): # 
    """
    Extracts the time from a given line.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    float: The extracted time.
    
    Raises:
    IndexError: If no time is found in the line.
    
    """
    try:
        return float(re.findall(r'(\d+.\d+)\s+Hours',line)[0])
    except:
        pass
def get_axial_offset(line:str):
    """
    Extracts the axial offset from a given line.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    float: The extracted axial offset.
    
    Raises:
    IndexError: If no axial offset is found in the line.
    
    """
    return float(re.findall(r'A-O\s+(\d+.\d+)',''.join(line))[0])
def get_boron_concentration(line:str):
    """
    Extracts the boron concentration from a given line.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    float: The extracted boron concentration.
    
    Raises:
    IndexError: If no boron concentration is found in the line.
    
    """
    return float(re.findall(r'BOR\s+(\d+.\d+)',''.join(line))[0])
def check_if_any_zero(numbers):
    for number in numbers:
        if number == 0:
            return True
    return False
def check_output_summary(line:str):
    """
    Checks if the the program has reached the Output summary and activates the output summary search.
    
    Parameters:
    line (str): A string containing the line of data.
    
    Returns:
    None
    """
    pattern = r"Output Summary"
    if re.search(pattern,line):
        Summary_search_activate()
        print("Output summary found")
    else:
        pass

#class cases containing the data for each case
class Case:
    """
    A class for storing the data of a case.
    
    Attributes:
    case_number (int): The case number.
    case_name (str): The case name.
    SA1 (XY_Nodal_Values): The fast absoption cross section nodal values.
    SA2 (XY_Nodal_Values): The thermal absorption cross section nodal values.
    NF1 (XY_Nodal_Values): The fast fission cross sectio nodal values.
    NF2 (XY_Nodal_Values): The thermal fission cross section nodal values.
    D1 (XY_Nodal_Values): The fast diffusion coefficient nodal values.
    D2 (XY_Nodal_Values): The thermal diffusion coefficient nodal values.
    KF1 (XY_Nodal_Values): The fast kappa*fission cross section nodal values.
    KF2 (XY_Nodal_Values): The thermal kappa*fission cross section nodal values.
    SS12 (XY_Nodal_Values): The fast to thermal scattering cross section nodal values.
    FLX1 (FLX_Nodal_Values): The fast neutron flux nodal values.
    FLX2 (FLX_Nodal_Values): The thermal neutron flux nodal values.
    XENON (XZ_Nodal_Values): The xenon concentration nodal values.
    IODINE (XZ_Nodal_Values): The iodine concentration nodal values.
    """
    def __init__(self,case_number:int):
        self.case_number = case_number
        self.case_name = "Case_"+str(case_number)
        self.SA1 = XY_Nodal_Values("SA1")
        self.SA2 = XY_Nodal_Values("SA2")
        self.NF1 = XY_Nodal_Values("NF1")
        self.NF2 = XY_Nodal_Values("NF2")
        self.D1 = XY_Nodal_Values("D1")
        self.D2 = XY_Nodal_Values("D2")
        self.KF1 = XY_Nodal_Values("KF1")
        self.KF2 = XY_Nodal_Values("KF2")
        self.SS12 = XY_Nodal_Values("SS12")
        self.FLX1 = FLX_Nodal_Values("FLX - Group 1")
        self.FLX2 = FLX_Nodal_Values("FLX - Group 2")
        self.XENON = XZ_Nodal_Values("XENON-135")
        self.IODINE = XZ_Nodal_Values("IODINE-135")
        self.Power = 0 # relative power in MW
        self.Rel_Power = 0 # Relative power in % 
        self.Boron_conc = 0 # Boron concentration in ppm
        self.AO = [] # Axial offset
        self.Exposure = 0 # core average exposure in GWd/MT
        self.TAP = [] #Time After Perturbation
    def check_stationary_variables(self):
        numbers = [self.Rel_Power,self.Power,self.Boron_conc,self.Exposure]
        return check_if_any_zero(numbers)
        

    # write a method to get only the names of XY_Nodal_Values
    def get_XY_Nodal_Values_name_list(self):
        """
        Gets the names of all the XY_Nodal_Values in the case.

        Returns:
        list: A list of the names of all the XY_Nodal_Values in the case.
        
        """
        xy_nodal_values_names = [
            name for name in dir(self)
            if not name.startswith('__') and isinstance(getattr(self, name), XY_Nodal_Values)
        ]
        return xy_nodal_values_names
    def get_XZ_Nodal_Values_name_list(self):
        """
        Gets the names of all the XZ_Nodal_Values in the case.
        
        Returns:
        list: A list of the names of all the XZ_Nodal_Values in the case.
        
        """
        xz_nodal_values_names = [
            name for name in dir(self)
            if not name.startswith('__') and isinstance(getattr(self, name), XZ_Nodal_Values)
        ]
        return xz_nodal_values_names
    def get_XZ_variable_match_list(self):
        """
        Gets the variable match strings of all the XZ_Nodal_Values in the case.
        
        Returns:
        list: A list of the variable match strings of all the XZ_Nodal_Values in the case.
        
        """
        xz_variable_match_list = [
            getattr(self, name).variable_match for name in dir(self)
            if not name.startswith('__') and isinstance(getattr(self, name), XZ_Nodal_Values)
        ]
        return xz_variable_match_list
    def get_name_list(self):
        """
        Gets the names of all the Nodal_Values in the case.
        
        Returns:
        list: A list of the names of all the Nodal_Values in the case."""
        return self.get_XY_Nodal_Values_name_list()+self.get_XZ_Nodal_Values_name_list()


# All nodal values are 3 dimensional arrays and have the same dimensions. They differ in the way they are presented in the .OUT file
class Nodal_Values:
    """
    A class for storing nodal values.
    
    Attributes:
    variable_match (str): The variable match string of the nodal values.
    
    """
    def __init__(self,variable_match:str):
        self.variable_match = variable_match
        self.data = np.full((no_radial_nodes,no_radial_nodes,no_axial_nodes),np.nan)
        self.name = self.variable_match
    def get_x_indices(self,line:str):
        """
        Gets the x indices of a line.
        
        Parameters:
        line (str): A string containing the line of data.
        
        Returns:
        list: A list of the x indices of the line.
        
        """
        return [int(i) for i in re.findall(r'\d\d*',"".join(line))]
    def get_nodal_data(self,line:str):
        """
        Gets the nodal data from a line.
        
        Parameters:
        line (str): A string containing the line of data.
        
        Returns:
        list: A list of the nodal data from the line.
        
        """
        return [float(i) for i in re.findall(r'\d+\.\d+',"".join(line))]
    def get_scaling(self,line:str): # Gets the scaling of the nodal values. if the scaling is 1 no value is shown.
        """
        Gets the scaling of the nodal values.
        
        Parameters:
        line (str): A string containing the line of data.
        
        Returns:
        float: The scaling of the nodal values.
        
        """
        if "PRI" in line:
            try:
                return float(re.findall(r"REAL VALUE \* (.*) = EDIT VALUE",''.join(line))[0])
            except:
                return 1.0
        else:
            error("Line passed to get_scaling() did not contain \"PRI\".\n Line passed to get_scaling() was: \n"+line) 
# Cross sections are presented in the XY plane and pages through K values
class XY_Nodal_Values(Nodal_Values):
    """
    A class for storing XY nodal values ie. macroscopic Cross sections.
    
    Attributes:
    variable_match (str): The variable match string of the nodal values.
    
    """
    def __init__(self, variable_match: str):
        super().__init__(variable_match)
        self.re_variable_match = re.compile(r"PRI.MAC - Cross section "+variable_match)
    def get_y_index(self,line:str):
        """
         Gets the y index of a line.
         
         Parameters:
         line (str): A string containing the line of data.
         
         Returns:
         int: The y index of the line.
         
         """
        return int(re.findall(r'\d\d?',"".join(line))[0])
    def get_z_index(self,line:str):
        """
         Gets the z index of a line.
         
         Parameters:
         line (str): A string containing the line of data.
         
         Returns:
         int: The z index of the line.
         
         """
        return int(re.findall(r'K\s+=\s+(\d\d?)',"".join(line))[0])
    
# Neutron fluxes and fission product concentrations are shown in the XZ plane and pages in the Y direction.
class XZ_Nodal_Values(Nodal_Values):
    """
    A class for storing XZ nodal values ie., Xenon and Iodine concentration as well as the flux.
    
    Attributes:
    variable_match (str): The variable match string of the nodal values.
    
    """
    def __init__(self, variable_match: str):
        super().__init__(variable_match)
        self.re_variable_match = re.compile(r"PRI.STA \w+\s+-\s+NODAL 3D "+variable_match)
        self.name = " ".join(re.findall("[a-zA-Z]+", variable_match))
        #Initialise time array
        self.data_array = np.zeros(list_length*n_timesteps)
    def get_y_index(self,line:str):
        return int(re.findall(r'\( *(\d+)\)',''.join(line))[0])
    def get_z_index(self,line:str):
        return int(re.findall(r'\d\d?',"".join(line))[0])
    
    
class FLX_Nodal_Values(XZ_Nodal_Values):
    """
    A class for storing FLX nodal values.
    
    Attributes:
    variable_match (str): The variable match string of the nodal values.
    
    """
    def __init__(self, variable_match: str):
        super().__init__(variable_match)
        self.group = int(re.findall(r'(\d+)',variable_match)[0])
        # the name should be FLX1 or FLX2
        self.name = "FLX"+str(self.group)

    
if __name__ == "__main__":
    # Get the file path
    myfile = get_out_file_path()
    #Initialise global variables
    initialise_global_variables(myfile)
    file_name = os.path.splitext(os.path.basename(myfile))[0]
    #Initialise dictionary of cases
    mydict = {}
    for case in cases:
        instance = Case(case)
        print("Initialising "+instance.case_name)
        for variable in [a for a in dir(instance) if not a.startswith('__')]:
            if isinstance(eval("instance."+variable),XY_Nodal_Values) or isinstance(eval("instance."+variable),XZ_Nodal_Values):
                mydict["Case_"+str(case)+"_"+variable] = eval("instance."+variable+".data")
            else:
                mydict["Case_"+str(case)+"_"+variable] = eval("instance."+variable)
        print(instance.case_name+" initialised")
    #Loop through lines and extract data
    Case_search_activate()
    line_number = 0
    for line in lines:
        line_number += 1
        check_output_summary(line)
        if  read_test_FLX1 == FALSE:   
            if ((mydict["Case_2_FLX1"][:,17,0] != test_FLX1).any()):
                    Fast_flux_changed = 1
            if Fast_flux_changed:
                print("Current case number = " + str(current_case.case_number) +"\n Current variable name = " + current_variable.name +  "\n Line number = " + str(line_number))
                error("Fast flux for case 2 has changed")
            #else:
                #print("Current case number = " + str(current_case.case_number) +"\n Current variable name = " + current_variable.name + "\n Line number = " + str(line_number))
                #print("Fast flux for case 2 hasnt changed")
        if search_case: #When search_case is true the search for the case number is active
            if case_line_identifier(line):
                current_case_number = get_case_number(line)
                current_case = Case(current_case_number)
                current_step_number = get_step_number(line)
                current_time = get_time(line)
                if (re.search(r"PRI.MAC\s+-\s+Cross section",lines[line_number+1]) or re.search(r"PRI.STA \w+\s+-\s+NODAL 3D",lines[line_number+1])):
                    Variable_search_activate()   
                continue
        elif (re.search(r"PRI.MAC\s+-\s+Cross section",line) and search_case == FALSE): #This search should get the cross sections in XY plane
            for variable in current_case.get_XY_Nodal_Values_name_list():
                if variable_line_identifier(line,variable): # Figure out which variable is being loaded
                    current_variable = eval("current_case."+variable)
                    current_scaling = current_variable.get_scaling(line)
                    current_z_index = current_variable.get_z_index(line)
                    Indices_search_activate() #When search_indices is true the x indices will be loaded.
                    break # End loop through variables when variable is found
            continue
        elif (re.search(r"PRI.STA \w+\s+-\s+NODAL 3D",line) and search_case == FALSE):
            for variable_match in current_case.get_XZ_variable_match_list(): # Figure out which variable is being loaded
                if variable_line_identifier(line,variable_match):
                    #Get the index of the match in the list of variable matches
                    idx = current_case.get_XZ_variable_match_list().index(variable_match)
                    current_variable = eval("current_case."+current_case.get_XZ_Nodal_Values_name_list()[idx])
                    current_scaling = current_variable.get_scaling(line)
                    Indices_search_activate()# When search_indices is true the x indices will be loaded.
                    search_y_index = TRUE # When search_y_index is true the y index will be loaded.
                    break # End loop through variables when variable is found
            continue      
        if search_indices:
            current_x_indices = current_variable.get_x_indices(line)
            current_x_indices = [item - 1 for item in current_x_indices]
            search_indices = FALSE
            search_data = TRUE # When search data is active the data will be loaded
            continue
        if search_data: 
            if isinstance(current_variable,XY_Nodal_Values):
                if line.strip() == "":
                    continue
                try:
                    if (current_y_index is not NONE and current_y_index == no_radial_nodes):
                        current_y_index = 1
                except NameError:
                    pass
                current_y_index = current_variable.get_y_index(line)
                current_line_data = current_variable.get_nodal_data(line)
                current_length = len(current_line_data)
                if current_x_indices[0] == 0:
                    current_variable.data[current_x_indices[-current_length:],current_y_index-1,current_z_index] = [item / current_scaling for item in current_line_data]
                    if current_y_index == no_radial_nodes: # This the last line for this half of the symmetry in this z index
                        search_data = FALSE
                        search_indices = TRUE
                        continue
                else:
                    current_variable.data[current_x_indices[0:current_length],current_y_index-1,current_z_index] =  [item / current_scaling for item in current_line_data]
                    if current_y_index == no_radial_nodes: # This is the last line of data for this z index
                        Case_search_activate() # When search_data is true the data will be loaded
                        mydict[current_case.case_name+"_"+current_variable.variable_match][:,:,current_z_index] = current_variable.data[:,:,current_z_index]
                        if (current_variable.variable_match == "SS12" and current_z_index == no_axial_nodes-1): # This is the last line of data for the variable
                            Case_search_activate() # When search_case is true the search for the case number is active
                        continue    
            if isinstance(current_variable,XZ_Nodal_Values):
                if search_y_index:
                    current_y_index = current_variable.get_y_index(line)
                    search_y_index = FALSE
                    search_data = TRUE
                    continue
                if search_data:
                    current_z_index = current_variable.get_z_index(line)
                    current_line_data = current_variable.get_nodal_data(line)
                    current_length = len(current_line_data)                    
                    if current_x_indices[0] == 0:
                        if current_y_index % 2 == 1:
                            current_variable.data[current_x_indices[-current_length:],current_y_index-1,current_z_index] = [item / current_scaling for item in current_line_data]
                            if current_z_index == 0: # This the last line for this half of the symmetry in this z index
                                search_data = TRUE
                                search_y_index = TRUE
                                search_x_indices = TRUE
                                continue
                        else:
                            current_variable.data[current_x_indices[-current_length:],current_y_index-1,current_z_index] = [item / current_scaling for item in current_line_data]
                            if current_z_index == 0: # This the last line for this half of the symmetry in this z index
                                search_data = FALSE
                                search_y_index = TRUE
                                search_x_indices = TRUE
                                continue
                    else:
                        if current_y_index % 2 == 1:
                            current_variable.data[current_x_indices[0:current_length],current_y_index-1,current_z_index] =  [item / current_scaling for item in current_line_data]
                            if current_z_index == 0: # This is the last line of data for this z index
                                if current_y_index == no_radial_nodes: # This is the last line of data for the variable
                                    mydict[current_case.case_name+"_"+current_variable.name] = current_variable.data
                                    Case_search_activate()
                                    search_x_indices = FALSE
                                    search_y_index = FALSE
                                    continue
                                search_data = TRUE
                                search_y_index = TRUE
                                serach_x_indices = TRUE
                                continue
                        else:
                            current_variable.data[current_x_indices[0:current_length],current_y_index-1,current_z_index] =  [item / current_scaling for item in current_line_data]
                            if current_z_index == 0: # This is the last line of data for this z index
                                if current_y_index == no_radial_nodes: # This is the last line of data for the variable
                                    mydict[current_case.case_name+"_"+current_variable.name] = current_variable.data
                                    if (read_test_FLX1 and current_case.case_number == 2):
                                        test_FLX1 = mydict["Case_2_FLX1"][:,17,0]
                                        read_test_FLX1 = FALSE                                                                                                  
                                    Case_search_activate()
                                    search_x_indices = FALSE
                                    search_y_index = FALSE
                                    continue
                                search_data = FALSE
                                search_y_index = TRUE
                                serach_x_indices = TRUE
                                continue
        elif search_summary:
                if not current_case.check_stationary_variables(): # Checks if all the stationary values have been gathered
                    mydict[current_case.case_name+"_"+"Exposure"] = current_case.Exposure
                    mydict[current_case.case_name+"_"+"Boron_conc"] = current_case.Boron_conc
                    mydict[current_case.case_name+"_"+"AO"] = current_case.AO
                    mydict[current_case.case_name+"_"+"Rel_Power"] = current_case.Rel_Power
                    mydict[current_case.case_name+"_"+"Power"] = current_case.Power
                    mydict[current_case.case_name+"_"+"Full_Power"] = current_case.Power/(current_case.Rel_Power*0.01)
                    print("All Stationary variables written for "+ str(current_case.case_name))
                    Case_search_activate()
                    continue
                try:
                    current_case.Exposure = get_average_exposure(line)
                except IndexError:
                    pass
                try:
                    current_case.Boron_conc = get_boron_concentration(line)
                except IndexError:
                    pass
                try:
                    current_case.AO = get_axial_offset(line)
                except IndexError:
                    pass
                try:
                    current_case.Power = get_power_level(line)
                except IndexError:
                    pass
                try:
                    current_case.Rel_Power = get_relative_power_level(line)
                except IndexError:
                    pass
    mydict["time_array"] = time_array
    # Save data to .mat file
    print("Saving data to data.mat")    
    sio.savemat("input/"+file_name+".mat",mydict)
    print("code completed")
        