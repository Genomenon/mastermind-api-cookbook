import requests, urllib, sys

URL = "https://mastermind.genomenon.com/api/v2/"
# Find your API token by logging in, visiting https://mastermind.genomenon.com/api, and clicking the link that says "Click here to fetch your API token".
API_TOKEN = "INSERT API TOKEN HERE"

def json_or_print_error(response):
    if response.status_code == requests.codes.ok:
        return response.json()
    else:
        print(response.status_code)
        print(response.text)
        sys.exit(0)

def encode(str):
    if sys.version_info[0] < 3:
        return urllib.quote_plus(str)
    else:
        return urllib.parse.quote_plus(str)

def main():

    variant = raw_input("Type genomic variant (e.g. NC_000012.11:g.57489193T>C): ")

    response = requests.get(url=URL+"suggestions", params={'api_token': API_TOKEN, 'variant': variant})

    data = json_or_print_error(response)
    encoded_match = encode(data[0]['matched'])

    print("Matched: " + data[0]['matched'])
    print("Mastermind link: https://mastermind.genomenon.com/detail?disease=all%20diseases&mutation=" + encoded_match)
    print("Maps to " + str(len(data)) + " effects:")

    for match in data:
        canonical = match['canonical']
        count_response = requests.get(url=URL+"counts", params={'api_token': API_TOKEN, 'variant': canonical})
        count_data = json_or_print_error(count_response)
        print("\t" + canonical + " - " + str(count_data['article_count']) + " articles, direct link: " + count_data['url'])

if __name__ == "__main__":
    main()
