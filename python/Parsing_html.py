# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 15:06:48 2014

@author: Zhen
"""

from bs4 import BeautifulSoup
import requests
import sys   

state = sys.argv[1]
state = "CT"
outfile = "%s_nursing_homes.csv"%(state) #Name of output file
outf = open(outfile, "w")
outf.write('Name,Number of Beds,Number of Residents,City,State,Contact\n')
    

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

for n in url_n:
    page = requests.get(n)
    soup = BeautifulSoup(page.content)
    name = soup.find_all('h1', {"align": "center"})
    
    for n in name:
        add = n.contents[0].split("-")
        hos = add[0]
        city = add[1].split(",")[0]
        state = add[1].split(",")[1]

    if soup.find("Number of Certified Beds:") == True:        
        info = soup.find_all('div', {"style": "float:left;padding-right:2.5em;"})
        for item in info:
            try:
                phone = item.contents[4].text
            except:
                phone = " "
            try:
                beds = item.contents[8].text
            except:
                beds = " "
            try:
                residents = item.contents[11].text
            except:
                residents = " "  
        outf.write('%s,%s,%s,%s,%s,%s\n' %(hos, 
                                       beds, 
                                       residents, 
                                       city,
                                       state,
                                       phone))
    #if  
        
#-----close file
outf.close()