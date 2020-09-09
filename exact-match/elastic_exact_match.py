from typing import List

from elasticsearch import Elasticsearch


def generate_must(fingerprint: str):
    bits = fingerprint.split(" ")
    json = []
    for bit in bits:
        json.append(
            {"term": {"fingerprint.keyword": {"value": "{}".format(bit)}}}
        )
    return json


def exact_match(
    es: Elasticsearch, fingerprint: str, index: str = "pubchem"
) -> List[str]:
    query = {"query": {"bool": {"must": generate_must(fingerprint)}}}

    res = es.search(body=query, index=index)
    hits = res["hits"]["hits"]
    if len(hits) > 0:
        result = []
        for hit in hits:
            result.append(hit["_source"]["smiles"])
        return result
    else:
        return []


if __name__ == "__main__":
    es = Elasticsearch(
        ["https://admin:admin@localhost:9200"],
        verify_certs=False,
        ssl_show_warn=False,
    )
    print(
        exact_match(
            es,
            "1629 1638 1646 1648 1660 1661 1670 1673 1695 1697 1698 1705 1709 1719 1743 1748 1749 1753 1771 1772 1790 1793 1806 1810 1820 1836 1838 1854 1856 1864 1869 1879 1880 1890 1900 1902 1909 1910 1914 1917 1928 1932 1950 1952 1963 1966 1971 1973 1983 1985 1994 2001 2004 2010 2018 2032 2043 2046 2049 2053 2056 2082 2095 2098 2106 2108 2110 2111 2130 2134",
        )
    )
