import requests

# Define the API endpoint URL provided.

import requests
import json
#users can only select from the following PTM models:
modeloptions = ["Phosphoserine_Phosphothreonine",
                "Phosphotyrosine",
                "N-linked_glycosylation",
                "O-linked_glycosylation",
                "Ubiquitination",
                "SUMOylation",
                "N6-acetyllysine",
                "Methylarginine",
                "Methyllysine",
                "Pyrrolidone_carboxylic_acid",
                "S-palmitoyl_cysteine",
                "Hydroxyproline",
                "Hydroxylysine"]

#http://www.musite.net:5000/musitedeep/Phosphoserine_Phosphothreonine;O-linked_glycosylation/MERSPAVCCQDPRAELVERVAAISVAHLEEAEEGPEPASNGVDPPPRARAASVIPGSASRPTPVRPSLSARKFSLQERPAGSCLEAQVGPYSTG
model=modeloptions[0]+";"+modeloptions[3] #for multiple models
url = ("https://api.musite.net/musitedeep/Phosphoserine_Phosphothreonine;O-linked_glycosylation/MVSKGEEDNMASLPATHELHIFGSINDVDFDMVGQGTGNPNEGYEELNLKSTKGDLQFSPWILVPHIGYGFHQYLPYPDGMSPFQAAMVDGSGYQVHRTMQFEDGASLTVNYRYTYEGSHIKGEAQVIGTGFPADGPVMTNTLTAADWCMSKMTYPNDKTIISTFKWSYITVNGKRYRSTARTTYFAKPMAANYLKNQPMYVFRKTELKHSM")
xurl = "http://www.musite.net:5000/musitedeep/"+model+"/"
myResponse = requests.get(url)
if(myResponse.ok):
       # In this Example, jData are prediction results from MusiteDeep predictor
       jData = json.loads(myResponse.content.decode('utf-8'))
       if "Error" in jData.keys(): 
           print(jData["Error"])
       else:
           print(jData)
else:
    myResponse.raise_for_status()
