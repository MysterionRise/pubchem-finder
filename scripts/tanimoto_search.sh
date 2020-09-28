GET pubchem/_search
{
  "query": {
    "script_score": {
      "query": {
        "bool": {
          "should": [
            {
              "term": {
                "fingerprint": {
                  "value": 1
                }
              }
            },
            {
              "term": {
                "fingerprint": {
                  "value": 2
                }
              }
            },
            {
              "term": {
                "fingerprint": {
                  "value": 3
                }
              }
            },
            {
              "term": {
                "fingerprint": {
                  "value": 4
                }
              }
            }
          ]
        }
      },
      "script": {
        "source": "_score / (params.a + doc['fingerprint'].length - _score) ",
        "params": {
          "a": 4
        }
      }
    }
  }
}