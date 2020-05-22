import sys
import requests
import urllib

import settings

def api_request(endpoint, options, request_type="GET", json_request=True, tries=0):
    params = options.copy()
    params.update({'api_token': settings.API_TOKEN})

    # print("Querying API: ", endpoint, options)
    response = requests.request(request_type, url=settings.URL+endpoint, params=params)

    if json_request:
        return json_or_print_error(response, endpoint, options, request_type, json_request, tries)
    else:
        return response

def json_or_print_error(response, endpoint, options, request_type, json_request, tries):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        sys.stdout.write('\n')
        print("ERROR ENCOUNTERED. ERROR CODE: " + str(response.status_code))
        if response.status_code != 500 and response.text:
            print("\t" + response.text)
        print("\tRESULTING FROM REQUEST: " + endpoint)
        print("\tWITH PARAMS: " + str(options))
        if response.status_code in [408, 500]:
            # Only retry GET requests
            if tries < 1 and request_type == "GET":
                print("Time out error, trying again.")
                return api_request(endpoint, options, request_type, json_request, tries+1)
            else:
                print("SKIPPING DATA FOR ABOVE REQUEST")
                return
        sys.exit(0)

def encode(str):
    if sys.version_info[0] < 3:
        return urllib.quote_plus(str)
    else:
        return urllib.parse.quote_plus(str)

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = 'M' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()
