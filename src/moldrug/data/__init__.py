from moldrug.data._get import get_data
import json

if __name__ == "__main__":

    print(json.dumps(get_data('6lu7'), indent=3))
