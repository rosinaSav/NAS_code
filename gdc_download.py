import generic as gen
import requests
import json

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
    "size": "200000"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

data = response.content.decode("utf-8")

json_data = json.loads(data)

print(json_data["data"]["pagination"])

ids = []
for hit in json_data["data"]["hits"][:3]:
    filename = hit["file_name"]
    cases = hit["cases"]
    id = hit["id"]
    ids.append(id)

outfile_path = "temp_data/json_data.txt"
outfile = open(outfile_path, "w")
id_list = {"ids": ids}
outfile.write(json.dumps(id_list))
outfile.close()

# str = "curl --remote-name --remote-header-name --request POST --header 'Content-Type: application/json' --data @temp_data/json_data.txt 'https://api.gdc.cancer.gov/legacy/data'"
# print(str.split(' '))
args = ['curl', '--remote-name', '--remote-header-name', '--request', 'POST', '--header', '"Content-Type: application/json"', '--data', '@{0}'.format(outfile_path), '"https://api.gdc.cancer.gov/legacy/data"']
" ".join(args)
gen.run_process(args)
