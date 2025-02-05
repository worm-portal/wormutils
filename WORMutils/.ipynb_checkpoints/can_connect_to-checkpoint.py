import requests

def can_connect_to(url) -> bool:  # `bool` is a type hint for the return type of the function
    try:
        response = requests.get(url)
    except requests.ConnectionError:
        return False  # can NOT connect
    else:
        if response.status_code == 200:
            return True  # can connect
        else:
            print("Unanticipated error code when attempting to reach a URL ("+str(url)+"):", response.status_code)
            return False