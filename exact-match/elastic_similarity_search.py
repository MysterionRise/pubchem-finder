from typing import List
import json
from elasticsearch import Elasticsearch


def generate_clauses(bits: List[int]):
    json = []
    for bit in bits:
        json.append({
            "term": {
                "fingerprint": {
                    "value": "{}".format(bit)
                }
            }
        })
    return json


def tversky_match(es: Elasticsearch, fingerprint: str, alpha: float, beta: float, index: str = "pubchem") -> List[str]:
    bits = fingerprint.split(' ')
    query = {
        'query': {
            'script_score': {
                'query': {
                    'bool': {
                        'should': generate_clauses(bits)
                    }
                },
                'script': {
                    'source': "_score / ((params.a - _score) * params.alpha + (doc['fingerprint'].length - _score) * params.beta)",
                    'params': {
                        'a': len(bits),
                        'alpha': alpha,
                        'beta': beta
                    }
                }
            }
        }
    }

    print(json.dumps(query))
    res = es.search(
        body=query,
        index=index
    )
    hits = res['hits']['hits']
    if len(hits) > 0:
        result = []
        for hit in hits:
            result.append((hit['_source']['smiles'], hit['_score']))
        return result
    else:
        return []


def tanimoto_match(es: Elasticsearch, fingerprint: str, index: str = "pubchem") -> List[str]:
    bits = fingerprint.split(' ')
    query = {
        'query': {
            'script_score': {
                'query': {
                    'bool': {
                        'should': generate_clauses(bits)
                    }
                },
                'script': {
                    'source': "_score / (params.a + doc['fingerprint'].length - _score) ",
                    'params': {
                        'a': len(bits)
                    }
                }
            }
        }
    }

    print(json.dumps(query))
    res = es.search(
        body=query,
        index=index
    )
    hits = res['hits']['hits']
    if len(hits) > 0:
        result = []
        for hit in hits:
            result.append((hit['_source']['smiles'], hit['_score']))
        return result
    else:
        return []
