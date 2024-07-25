import json
import sys
import csv
from urllib.error import HTTPError
from urllib.request import urlopen


def query_api(query):
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/all/protein/UniProt/{query}/?page_size=200&extra_fields=hierarchy,short_name"

    try:
        with urlopen(url) as res:
            data = json.loads(res.read().decode("utf-8"))
        return data

    except HTTPError as e:
        print(f"Error querying API for {query}: {e}")
        return None

    except Exception as e:
        print(f"Unexpected error querying API for {query}: {e}")
        return None


def extract_features(data):
    features_data = []

    for entry in data["results"]:
        meta = entry["metadata"]
        protein = entry["proteins"][0]
        signatures = ",".join([",".join(db.keys()) for db in meta["member_databases"].values()]) if meta[
            "member_databases"] else "-"
        go_terms = ",".join([t["identifier"] for t in meta["go_terms"]]) if meta["go_terms"] else "-"
        locations = [f"{f['start']}..{f['end']}" for protein in data["results"] for l in
                     protein.get("entry_protein_locations", []) for f in l["fragments"]]

        features_data.append([
            meta["accession"],
            meta["name"] or "-",
            meta["source_database"],
            meta["type"],
            meta["integrated"] or "-",
            signatures,
            go_terms,
            protein["accession"].upper(),
            str(protein["protein_length"]),
            ",".join(locations)
        ])

    return features_data


def query_and_save(query_csv, output_csv):
    with open(query_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            query = row[0]
            data = query_api(query)
            features_data = extract_features(data)

            with open(output_csv, mode='a', newline='') as output_file:
                writer = csv.writer(output_file)
                for feature_row in features_data:
                    writer.writerow(feature_row)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_csv output_csv")
        sys.exit(1)

    query_csv = sys.argv[1]
    output_csv = sys.argv[2]
    query_and_save(query_csv, output_csv)
