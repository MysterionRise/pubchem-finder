from typing import List

from elasticsearch import Elasticsearch


def generate_must(fingerprint: str):
    bits = fingerprint.split(' ')
    json = []
    for bit in bits:
        json.append({
            "term": {
                "fingerprint.keyword": {
                    "value": "{}".format(bit)
                }
            }
        })
    return json


def exact_match(es: Elasticsearch, fingerprint: str, index: str = "pubchem") -> List[str]:
    query = {
        'query': {
            'bool': {
                'must': generate_must(fingerprint)
            }
        }
    }

    res = es.search(
        body=query,
        index=index
    )
    hits = res['hits']['hits']
    if len(hits) > 0:
        result = []
        for hit in hits:
            result.append(hit['_source']['smiles'])
        return result
    else:
        return []
