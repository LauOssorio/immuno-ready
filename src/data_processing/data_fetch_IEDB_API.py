import requests
import json
import pandas as pd
import time
from io import StringIO
from datetime import datetime
from src.config import *

base_uri= IEDB_API_BASE_URL


def print_curl_cmd(req):
    """
    Print a cURL command equivalent to the given request for debugging.
    """
    url = req.url
    print("curl -X 'GET' '" + url + "'")

def retrieve_IEDB_api_data():
    """
    Retrieve data from the IEDB API using pagination for requests larger than 10000
    observations and return it as a DataFrame.
    """

    search_params={'host_organism_iri_search': 'cs.{"NCBITaxon:9606"}',
                    'order': 'structure_id', 'offset':0,
                    'qualitative_measure': 'not.eq.Negative'
                    }

    table_name='tcell_search'

    full_url=base_uri + '/' + table_name

    result = requests.get(full_url, params=search_params)

    print_curl_cmd(result)

    df = pd.json_normalize(result.json())

    while len(result.json()) > 0:
        time.sleep(2)
        search_params['offset'] += 10000
        print('offset: ' + str(search_params['offset']))
        result = requests.get(full_url, params=search_params)
        df = pd.concat([df, pd.json_normalize(result.json())], ignore_index=True)

    print('Done!')

    return df


def api_request_to_csv(filename_prefix = "IEDB_positive_peptides"):
    """
    Save retrieved IEDB data to a CSV file in the raw data folder.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.csv"
    retrieve_IEDB_api_data().to_csv(RAW_DATA_PATH + filename, index=False)


if __name__ == "__main__":
    api_request_to_csv()
