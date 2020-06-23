###########################################
# BLAT
# Author: Dong Liang
# Email: ldifer@gmail.com
# May 20, 2020
###########################################


from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from bs4 import BeautifulSoup
import bs4
import pandas as pd
import numpy as np
import re
import os


class BLAT(object):
    def __init__(self, base_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"):
        self.base_url = base_url


    def open(self, headless = True):
        options = webdriver.ChromeOptions()
        if headless:
            options.add_argument("headless")
        driver = webdriver.Chrome(executable_path='/usr/local/bin/chromedriver', chrome_options=options)
        driver.get(self.base_url)
        driver.implicitly_wait(100)
        return driver
        # driver.get("https://www.google.com")
        # driver.get_screenshot_as_file("~/Desktop/capture.png")


    def query(self, query, species = 'human', headless = False):
        # Open the weblink
        driver = self.open(headless = headless)

        # Select species
        popup = driver.find_element_by_name('org')
        popup.send_keys(species)

        # Submit the query
        search_box = driver.find_element_by_name('userSeq')
        search_box.send_keys(query)
        search_box.submit()

        # Parse the returned html page
        html = driver.page_source
        soup = BeautifulSoup(html, 'html.parser')
        content = soup.find('pre')

        blat_res = list()
        browser = list()
        details = list()
        
        for i, c in enumerate(content):
            # Extract text contents 
            if isinstance(c, bs4.element.NavigableString) and i != 0 and len(c.split()) != 0:
                blat_res.append(c.split())

            # Extract weblinks (href)
            if isinstance(c, bs4.element.Tag):
                if c.contents[0] == 'browser':
                    browser.append('https://genome.ucsc.edu/' + c.get('href'))
                elif c.contents[0] == 'details':
                    details.append('https://genome.ucsc.edu/' + c.get('href'))
        
        # Close the chromedriver
        driver.close()

        # Wrap up in a dataframe
        blat_res = np.array(blat_res)
        columns = ['BROWER', 'DETAILS', 'QUERY', 'SCORE', 'START', 'END', 'QSIZE', 'IDENTITY', 'CHROM', 'STRAND', 'START', 'END', 'SPAN']
        blat_res = pd.DataFrame(np.c_[browser, details, blat_res], columns = columns)

        return blat_res

    
        

import csv
class BLAT_summary(object):
    def __init__(self):
        # self.filename = filename
        self.rcd = []

    def file_handler(self, infile):
        file = pd.read_excel(infile)
        file = file.query("chr_name == 'chr8' and mapq > 30 and alignment_len > 1000")
        seq = file.seq 
        header = file.header
        print("-- -- -- -- --")
        print(infile + "..." + str(len(seq)))
        return seq, header

  
    
    def fetch_blat(self, seq, header):
        blat = BLAT()

        try:
            hg38 = blat.query(seq, 'human', True)
            self.rcd.append([str(header)] + ['BLAT hg38'] + list(hg38.iloc[0, :].values))
        except:
            self.rcd.append([str(header)] + ['BLAT hg38'] + [])

        try:
            mm10 = blat.query(seq, 'mouse', True)
            self.rcd.append([str(header)] + ['BLAT mouse'] + list(mm10.iloc[0, :].values))
        except:
            self.rcd.append([str(header)] + ['BLAT mouse'] + [])

        # self.rcd.append([str(header)] + ['BLAT hg38'] + list(hg38.iloc[0, :].values))
        # self.rcd.append([str(header)] + ['BLAT mouse'] + list(mm10.iloc[0, :].values))
        print(header, 'done')

        # with open('output.csv', mode = 'a+') as csv_file:
        #     writer = csv.writer(csv_file)
        #     # writer.writerow(['Header', 'BLAT species'] + ['BROWER', 'DETAILS', 'QUERY', 'SCORE', 'START', 'END', 'QSIZE', 'IDENTITY', 'CHROM', 'STRAND', 'START', 'END', 'SPAN'])
        #     writer.writerow([str(header)] + ['BLAT hg38'] + list(hg38.iloc[0, :].values))
        #     writer.writerow([str(header)] + ['BLAT mouse'] + list(mm10.iloc[0, :].values))
        # print('Write to output.csv')
        # return np.c_[mm10.iloc[0, :][2::], hg38.iloc[0, :][2::]]

    def __call__(self, infile):
        seq, header = self.file_handler(infile)
        output_filename = infile.split('.')[0]+ '_BLAT.csv'
        with open("./BLATb/" + output_filename, mode = 'w+') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Header', 'BLAT species'] + ['BROWER', 'DETAILS', 'QUERY', 'SCORE', 'START', 'END', 'QSIZE', 'IDENTITY', 'CHROM', 'STRAND', 'START', 'END', 'SPAN'])
            for s, h in zip(seq, header):
                self.fetch_blat(s, h)
            writer.writerows(self.rcd)

        print(f'Written to {output_filename}')
        


###########################################
# Hit
# Author: Dong Liang
# Email: ldifer@gmail.com
# May 23, 2020
###########################################


def hit(infile):
    rcd = pd.read_csv(infile)  

    # Proprocessing
    rcd.rename({'BLAT species': 'BLAT_species'}, axis='columns', inplace = True)
    rcd.IDENTITY = rcd.IDENTITY.apply(lambda x: float(str(x).strip('%')))
    rcd.drop(columns = ['BROWER', 'DETAILS'], inplace = True)
    
    # Condition
    condition1 = "SCORE > 0 and IDENTITY > 99 and END < QSIZE and END - START > 20" # hg38
    condition2 = "SCORE > 0 and IDENTITY > 99 and START != 1 and END - START > 20"   # mm10
    if any(rcd.query(condition1)['BLAT_species'] != 'BLAT hg38'):
        print(rcd.query(condition1))
    elif any(rcd.query(condition2)['BLAT_species'] != 'BLAT mouse'):
        print(rcd.query(condition2))
    

def hit2(infile):
    rcd = pd.read_csv(infile)  

    # Proprocessing
    rcd.rename({'BLAT species': 'BLAT_species'}, axis='columns', inplace = True)
    rcd.IDENTITY = rcd.IDENTITY.apply(lambda x: float(str(x).strip('%')))
    rcd.drop(columns = ['BROWER', 'DETAILS'], inplace = True)
    
    # Header
    header = rcd.Header.unique()
    both = list()
    pure = list()
    rcd['Type'] = ''
    for h in header:
        obs = rcd.loc[rcd.Header == h]
        if all(obs.IDENTITY > 99) and all(obs.SCORE > 100):
            rcd.loc[rcd.Header == h, 'Type'] = 'Both'
            
        elif any(abs(obs.END - obs.START - obs.QSIZE) + 1 < 10) :
            # print('Condition2: pure')
            rcd.loc[rcd.Header == h, 'Type'] = 'Pure'
            # pure.append(obs)            
        else:
            rcd.loc[rcd.Header == h, 'Type'] = 'Overlapped'
        
    
    print('#################### Condition1: both #####################')
    print(rcd.loc[rcd.Type == 'Both'])
    print('#################### Condition2: pure #####################')
    print(rcd.query("Type == 'Pure' and IDENTITY > 99")) # [rcd.Type == 'Pure' and rcd.IDENTITY > 99, :])
    print('#################### Condition2: Overlapped #####################')
    print(rcd.loc[rcd.Type == 'Overlapped'])


