### web.py ###
# Author: Marc Zepeda
# Date: 2024-08-17

# Import packages
import requests
import os
from selenium import webdriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select, WebDriverWait
from selenium.common.exceptions import StaleElementReferenceException, NoSuchElementException, UnexpectedAlertPresentException
import time
import io


# HTML Methods
''' get_html: Obtains html file from the internet
        url: website url for html file
        file: saved file name
        dir: saved file directory
    Dependencies: os,requests
'''
def get_html(url: str(),file=None,dir=None):
    
    try: # Send a GET request to the URL
        response = requests.get(url)
        if response.status_code == 200: # Check if the request was successful
            if file is not None and dir is not None: # Save html file if file and dir are specified
                with open(os.path.join(dir,file), 'w', encoding='utf-8') as html: # Write the HTML content to a file
                    html.write(response.text)
                print(f"HTML content successfully saved to {file}")
        else: print(f"Failed to retrieve the webpage. Status code: {response.status_code}")
        return response.text
    except Exception as e: print(f"An error occurred: {e}")

''' search: Returns & print blocks of text that contains substring
        url: website url
        sub: substring
        within: blocks of text that contain substring
    Dependencies: get_html()
'''
def search(url: str(), sub: str(), within='\n'):
    texts=[block for block in get_html(url=url).split(within) if sub in block] # Isolate blocks of text that contain substring
    for text in texts: print(f'"{sub}" found between "{within}": {text}') # Print texts
    return texts

# Selenium Methods
''' handle_dropdowns: Recursively handle the next dropdown in the list
        dc: dictionary
        prev_option: previous dropdown menu selection
        driver: Selenium WebDriver (Chrome)
        dropdown_ids: list of dropdown menu ids (strings)
        current_index: enables iteration through dropdown_ids
    Dependencies: selenium,time
'''
def handle_dropdowns(dc: dict(), prev_option: str(), driver, dropdown_ids, current_index=1):
    
    if current_index >= len(dropdown_ids): # Base case: If we've handled all dropdowns, return
        return

    try: # Enables the handling of exceptions
        dropdown = Select(WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, dropdown_ids[current_index])))) # Locate the current dropdown by its ID
        options = [option.text for option in dropdown.options] # Get all the options in the dropdown
        dc[prev_option] = {key: None for key in options} # Dictionary of current dropdown options
        
        for j,option in enumerate(options): # Iterate through all options in the current dropdown
            print(f'Previous option: {prev_option}\nCurrent Index: \n{current_index}')
            print(f'j: {j}\nOption: {option}')
            
            dropdown = Select(WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, dropdown_ids[current_index])))) # Reinitialize the Select object due to potential DOM changes
            dropdown.deselect_all() # Deselect previously selected options
            dropdown.select_by_index(j) # Select the current option by index
            
            time.sleep(1) # Wait until the next dropdown is loaded/updated before proceeding
            handle_dropdowns(dc[prev_option], option, driver, dropdown_ids, current_index + 1) # Recursively handle the next dropdown in the list
            
    except StaleElementReferenceException: # Handle the stale element by retrying the operation
        handle_dropdowns(dc, prev_option, driver, dropdown_ids, current_index)
    
    except NoSuchElementException: # Skip no such element
        print(f'NoSuchElementException:\nPrevious option: {prev_option}\nCurrent Index: \n{current_index}\nj: {j}\nOption: {option}')
    
    except UnexpectedAlertPresentException: # Handle the alert
        alert = driver.switch_to.alert
        alert_text = alert.text  # You can get the alert text if needed
        print(f'Alert detected with message: {alert_text}')
        alert.accept() # Accept the alert to dismiss it
    
    except NotImplementedError: # Handle not implemented error
        print(f"Not Implemented Error:\nPrevious option: {prev_option}\nCurrent Index: \n{current_index}\nj: {j}\nOption: {option}")

''' get_dropdown_options: Returns dictionary of dropdown menu options that depend on each other.
        url: website url
        dropdown_ids: list of dropdown menu ids (strings)
        dir: directory to save checkpoints (optional)
        file: file name to save checkpoints (optional)
        dc_checkpoint: dictionary from last save checkpoint (optional)
        dc_checkpoint: index of last save checkpoint (optional)
    Dependencies: selenium,time,tidy,io,handle_dropdowns()
'''
def get_dropdown_options(url: str(), dropdown_ids: list(), dir=None, file=None, dc_checkpoint={},dc_checkpoint_i=0):
    
    driver = webdriver.Chrome() # Set up the Selenium WebDriver (Chrome in this case)
    
    try:
        driver.get(url) # Navigate to the website
        time.sleep(3) # Wait for the page to load fully

        dropdown = Select(driver.find_element(By.ID, dropdown_ids[0])) # Locate the current dropdown by its ID
        options = [option.text for option in dropdown.options] # Get all the options in the dropdown
        if not dc_checkpoint: dc = {key: None for key in options} # Initialize dictionary of dropdown options (start)
        else: dc=dc_checkpoint # Initialize dictionary of dropdown options (checkpoint provided)

        for i,option in enumerate(options): # Iterate through all options in the current dropdown
            if i >= dc_checkpoint_i: # Skips previous dropdown options already stored in dc_checkpoint
                
                time.sleep(1) # Wait until the current dropdown is loaded/updated before proceeding
                dropdown = Select(WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, dropdown_ids[0])))) # Reinitialize the Select object due to potential DOM changes
                dropdown.deselect_all() # Deselect previously selected options
                dropdown.select_by_index(i) # Select the current option by index
                
                time.sleep(1) # Wait until the next dropdown is loaded/updated before proceeding
                handle_dropdowns(dc, option, driver, dropdown_ids, current_index=1) # Recursively handle the next dropdown in the list
                
                print(f'Completed loop {i} with "{option}":\n{dc}') # Status update
                if file is not None and dir is not None: io.save(dir=dir, file=f'{file}_{i}.csv', obj=t.dc_to_ls(dc)) # Save checkpoint
            

    finally: driver.quit() # Close the browser window

    return dc