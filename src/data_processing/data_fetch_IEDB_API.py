import pandas as pd
import requests
import pandas as pd

def fetch_epitopes_api(filters=None, select=None, order=None, batch_size=10000, table_name='tcell_search'):
    """
    Fetches all epitope records from the IEDB IQ API using POST requests with pagination.

    Args:
        filters (dict): Dictionary of filters for the query (e.g., {"hla_class": {"EQ": "HLA-I"}}).
        select (list): List of fields to select.
        order (list): List of fields to order by (required for pagination).
        batch_size (int): Number of records to fetch per request (max 10000).

    Returns:
        pd.DataFrame: DataFrame containing all retrieved records.
    """
    base_url = "https://query-api.iedb.org"
    table_name='epitope_search'
    url=base_url + '/' + table_name
    all_data = []
    offset = 0

    if filters is None:
        filters = {}
    if select is None:
        select = []
    if order is None:
        raise ValueError("You must specify an 'order' field for pagination.")

    while True:
        payload = {
            "filters": filters,
            "select": select,
            "order": order,
            "limit": batch_size,
            "offset": offset
        }

        response = requests.post(url, json=payload)
        response.raise_for_status()

        batch_data = response.json()

        if not batch_data:
            break  # No more data to fetch

        all_data.extend(batch_data)

        if len(batch_data) < batch_size:
            break  # Last batch received

        offset += batch_size

    return pd.DataFrame(all_data)

# Example usage:
if __name__ == "__main__":
    filters = {
        "hla_class": {"EQ": "HLA-I"}  # Example filter
    }
    select = ["structure_id", "linear_sequence", "hla_class"]
    order = ["structure_id"]

    df = fetch_epitopes_api(filters=filters, select=select, order=order)
    print(f"Retrieved {len(df)} records")
    # You can save it if you want:
    df.to_csv("data/raw/iedb_epitopes.csv", index=False)
