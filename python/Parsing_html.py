# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 15:06:48 2014

@author: Zhen
"""

from bs4 import BeautifulSoup
import requests
import sys
import time
import os   


if len(sys.argv) != 2:
    print "Please specify a state using abbreviation"
    sys.exit()


state = sys.argv[1] 
if len(state) != 2:
    print "Please use abbreviation of a state"
    print "Usage:"
    print "      python Parsing_html.py STATE_ABBREVIATION"
    sys.exit()
    
#state = "CT"
outfile = "%s_nursing_homes.csv"%(state) #Name of output file
outf = open(outfile, "w")
outf.write('Name,Number of Beds,Number of Residents,City,State,Contact\n')
outf.close()
print "File saved to %s"%os.getcwd()
print "---------------------------------------"
#----Cook the soup
url = "http://www.hospital-data.com/dirs/%s-Nursing-Homes.html"%(state)
page = requests.get(url)
soup = BeautifulSoup(page.content)
    
#----Get all the urls to all cities in a state that have nursing homes
city = []    
for link in soup.find_all('a', href = True):
    if "city/" in link['href']:
        city.append(link['href'])

nursing_homes = []
for c in city:
    url_c = "http://www.hospital-data.com/dirs/" + c
    page_c = requests.get(url_c)
    soup_c = BeautifulSoup(page_c.text)
    for link in soup_c.find_all('a', href = True):
        if "hospitals" in link['href']:
            nursing_homes.append(link['href'])

url_n = [] 
for n in nursing_homes:            
    url_n.append("http://www.hospital-data.com/" + n.lstrip("../../"))

outf = open(outfile, "a")
for n in url_n:
    page = requests.get(n)
    #page = requests.get("http://www.hospital-data.com/hospitals/DINAN-MEMORIAL-CENTER-TEDESCO-ANNEX-BRIDGEPORT.html")
    soup = BeautifulSoup(page.content)
    name = soup.find('h1', {"align": "center"})
    print "Retrieve information for %s"%name.text
    add = name.contents[0].split(" - ")
    hos = add[0].replace(',', '')
    city = add[1].split(",")[0]

    if "Number of Certified Beds:" in soup.text:        
        info = soup.find('div', {"style": "float:left;padding-right:2.5em;"})
        try:
            phone = info.contents[4].text
        except:
            phone = " "
        try:
            beds = info.contents[8].text
        except:
            beds = " "
        try:
            residents = info.contents[11].text
        except:
            residents = " "      
    elif "Number Of Beds" in soup.text:
        info = soup.find('div', {"style": "float:left;padding-right:2.5em;"})
        try:        
            phone = info.contents[4].text
        except:
            phone = " "
        info = soup.find('div', {"class": "fdivpos"})
        info = info.find_next('div', {"class": "fdivpos"})
        try:
            beds = info.contents[5].text
        except:
            beds = " "
        residents = " "
    
    
    outf.write('%s,%s,%s,%s,%s,%s\n' %(hos, 
                                       beds, 
                                       residents, 
                                       city,
                                       state,
                                       phone))
    print "Finished retrieving information for %s"%name.text
    print "---------------------------------------"
    time.sleep(5)

outf.close() #------close file
print "All information successfully retrieved for nursing_homes in %s"%state
    
        

