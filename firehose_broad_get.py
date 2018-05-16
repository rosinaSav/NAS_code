'''
Authors: Liam Abrahams.
Small wrapper to use firebrowse to dump json data
'''

import generic as gen
import firebrowse
import json

def main():

    description = "Download disease json."
    arguments = ["json_file"]
    args = gen.parse_arguments(description, arguments, flags = [])
    json_file = args.json_file

    # get the information from firehose
    data = firebrowse.Archives().StandardData(page_size=2000);
    data = json.loads(data)

    # write to file
    with open(json_file, "w") as outfile:
        json.dump(data, outfile)

if __name__ == "__main__":
    main()
