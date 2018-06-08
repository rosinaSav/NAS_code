import generic as gen
import requests
import json

fields = [
    "file_name",
    "cases.submitter_id",
    "cases.samples.sample_type",
    "cases.disease_type",
    "cases.project.project_id"
    ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/legacy/files"

# This set of filters is nested under an 'and' operator.
filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.program.name",
            "value": ["TCGA"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_category",
            "value": ["Simple nucleotide variation"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.experimental_strategy",
            "value": ["DNA-Seq"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_type",
            "value": ["Simple somatic mutation"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.platform",
            "value": ["Illumina HiSeq"]
            }
        },
        {
        "op": "=",
        "content":{
            "field": "files.access",
            "value": ["open"]
            }
        }
    ]
}

# A POST is used, so the filter parameters can be passed directly as a Dict object.
params = {
    "filters": filters,
    "fields": fields,
    "format": "json",
    "size": "2000000"
    }


def run_data_query(mutation_data_file):

    print("Querying GDC api...")
    # The parameters are passed to 'json' rather than 'params' in this case
    response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

    data = response.content.decode("utf-8")

    json_data = json.loads(data)

    ids = []
    for hit in json_data["data"]["hits"]:
        filename = hit["file_name"]
        cases = hit["cases"]
        id = hit["id"]
        ids.append(id)

    outfile = open(mutation_data_file, "w")
    id_list = {"ids": ids}
    outfile.write(json.dumps(id_list))
    outfile.close()


def main():

    description = "Get the information to download files from the GDC."
    args = gen.parse_arguments(description,  ["mutation_data_file"], flags = [])
    mutation_data_file = args.mutation_data_file

    run_data_query(mutation_data_file)

    args = ['curl', '--remote-name', '--remote-header-name', '--request', 'POST', '--header', '"Content-Type: application/json"', '--data', '@{0}'.format(mutation_data_file), '"https://api.gdc.cancer.gov/legacy/data"']
    print("\nCopy and run the command below:\n")
    print("{0}\n".format(" ".join(args)))


if __name__ == "__main__":
    main()
