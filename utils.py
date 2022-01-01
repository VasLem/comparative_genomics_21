import pandas as pd
def append_sheet_to_excel(name, sheet, df, index=False):
    with pd.ExcelWriter(name, engine='openpyxl', mode='a') as writer: 
        workBook = writer.book
        try:
            workBook.remove(workBook[sheet])
        except:
            print("Worksheet does not exist")
        finally:
            df.to_excel(writer, sheet_name=sheet,index=index)
            writer.save()

import requests
import re
def batch_cd_search(files):
    CDSEARCH_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"
    SEARCH_PARAMS = {
        "db": "cdd",
        "smode": "auto",
        "useid1": "true",
        "compbasedadj": "1",
        "filter": "true",
        "evalue": "3.0",
        "maxhit": "500",
        "dmode": "full",
        "tdata": "hits",
    }
    response = requests.post(CDSEARCH_URL, params=SEARCH_PARAMS, files=files)
    match = re.search(r"#cdsid\t(.+?)\n", response.text)
    if match:
        cdsid = match.group(1)
        return cdsid