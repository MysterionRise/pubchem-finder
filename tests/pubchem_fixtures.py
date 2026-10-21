"""Recorded-style PubChem payloads for the knowledge tools, and a respx router.

Shapes follow PUG REST / PUG View documentation and real responses; the live
tests (PUBCHEM_LIVE=1) check the same paths against the real service.
"""

import re

import httpx
import respx

REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"


def row(cid, title, smiles, formula="", mw="0"):
    return {
        "CID": cid,
        "Title": title,
        "SMILES": smiles,
        "ConnectivitySMILES": smiles,
        "MolecularFormula": formula,
        "MolecularWeight": mw,
    }


ASPIRIN = row(2244, "Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "C9H8O4", "180.16")
SALICYLIC = row(338, "Salicylic Acid", "C1=CC=C(C(=C1)C(=O)O)O", "C7H6O3", "138.12")
METHYL_SALICYLATE = row(
    4133, "Methyl Salicylate", "COC(=O)C1=CC=CC=C1O", "C8H8O3", "152.15"
)
BROKEN = row(99999, "Unparseable", "not-a-smiles")


def table(*rows):
    return {"PropertyTable": {"Properties": list(rows)}}


def fault(code, message, *details):
    return {"Fault": {"Code": code, "Message": message, "Details": list(details)}}


SYNONYMS = {
    "InformationList": {
        "Information": [
            {
                "CID": 2244,
                "Synonym": [
                    "aspirin",
                    "ACETYLSALICYLIC ACID",
                    "50-78-2",
                    "Aspirin",
                    "2-Acetoxybenzoic acid",
                    "50-78-3",
                    "Acetylsalicylate",
                ],
            }
        ]
    }
}

XREFS = {
    "InformationList": {
        "Information": [
            {
                "CID": 2244,
                "RegistryID": ["50-78-2", "CHEMBL25", "DB00945"],
                "PubMedID": [123, 456, 456],
                "PatentID": ["US-1234567-A", "EP-0000001-B1"],
            }
        ]
    }
}

DESCRIPTION = {
    "InformationList": {
        "Information": [
            {"CID": 2244, "Title": "Aspirin"},
            {
                "CID": 2244,
                "Description": "Aspirin is a member of the class of benzoic acids.",
                "DescriptionSourceName": "ChEBI",
                "DescriptionURL": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:15365",
            },
            {
                "CID": 2244,
                "Description": "x" * 3000,
                "DescriptionSourceName": "LiverTox",
            },
        ]
    }
}

ASSAY_COLUMNS = [
    "AID",
    "Panel Member ID",
    "SID",
    "CID",
    "Activity Outcome",
    "Target Accession",
    "Target GeneID",
    "Activity Value [uM]",
    "Activity Name",
    "Assay Name",
    "Bioassay Type",
    "PubMed ID",
    "RNAi",
]


def assay(aid, outcome, acc="", gene="", value="", name="", assay_name=""):
    return {
        "Cell": [
            aid,
            "",
            "845",
            "2244",
            outcome,
            acc,
            gene,
            value,
            name,
            assay_name,
            "",
            "",
            "",
        ]
    }


ASSAY_SUMMARY = {
    "Table": {
        "Columns": {"Column": ASSAY_COLUMNS},
        "Row": [
            assay(
                "100", "Active", "P23219", "5742", "1.67", "IC50", "COX-1 inhibition"
            ),
            assay("101", "Active", "P35354", "5743", "278", "IC50", "COX-2 inhibition"),
            assay(
                "102", "Active", "P23219", "5742", "0.5", "IC50", "COX-1 whole blood"
            ),
            assay(
                "102", "Active", "P23219", "5742", "0.7", "IC50", "COX-1 whole blood"
            ),
            assay("103", "Inactive"),
            assay("104", "Inconclusive"),
            assay("105", "Active", value="n/a"),
            assay("abc", "Active"),
        ],
    }
}

GHS_VIEW = {
    "Record": {
        "RecordType": "CID",
        "RecordNumber": 2244,
        "RecordTitle": "Aspirin",
        "Section": [
            {
                "TOCHeading": "Safety and Hazards",
                "Section": [
                    {
                        "TOCHeading": "Hazards Identification",
                        "Section": [
                            {
                                "TOCHeading": "GHS Classification",
                                "Information": [
                                    {
                                        "ReferenceNumber": 1,
                                        "Name": "Pictogram(s)",
                                        "Value": {
                                            "StringWithMarkup": [
                                                {
                                                    "String": "  ",
                                                    "Markup": [
                                                        {
                                                            "Start": 0,
                                                            "Length": 1,
                                                            "URL": "https://pubchem.ncbi.nlm.nih.gov/images/ghs/GHS07.svg",
                                                            "Type": "Icon",
                                                            "Extra": "Irritant",
                                                        }
                                                    ],
                                                }
                                            ]
                                        },
                                    },
                                    {
                                        "ReferenceNumber": 1,
                                        "Name": "Signal",
                                        "Value": {
                                            "StringWithMarkup": [{"String": "Warning"}]
                                        },
                                    },
                                    {
                                        "ReferenceNumber": 1,
                                        "Name": "GHS Hazard Statements",
                                        "Value": {
                                            "StringWithMarkup": [
                                                {
                                                    "String": "H302 (100%): Harmful if swallowed"
                                                },
                                                {
                                                    "String": "H315: Causes skin irritation"
                                                },
                                            ]
                                        },
                                    },
                                    {
                                        "ReferenceNumber": 1,
                                        "Name": "Precautionary Statement Codes",
                                        "Value": {
                                            "StringWithMarkup": [
                                                {
                                                    "String": "P264, P270, P301+P317, P330, P501 "
                                                    * 12
                                                }
                                            ]
                                        },
                                    },
                                    {
                                        "ReferenceNumber": 2,
                                        "Name": "Empty",
                                        "Value": {},
                                    },
                                    {
                                        "ReferenceNumber": 2,
                                        "Name": "Signal",
                                        "Value": {
                                            "StringWithMarkup": [{"String": "Danger"}]
                                        },
                                    },
                                ],
                            }
                        ],
                    }
                ],
            }
        ],
        "Reference": [
            {
                "ReferenceNumber": 1,
                "SourceName": "European Chemicals Agency (ECHA)",
                "URL": "https://echa.europa.eu/substance-information/-/substanceinfo/100.000.059",
                "LicenseNote": "Use is subject to ECHA's legal notice",
                "LicenseURL": "https://echa.europa.eu/web/guest/legal-notice",
            },
            {
                "ReferenceNumber": 2,
                "SourceName": "Some Depositor",
                "URL": "https://example.org",
            },
        ],
    }
}

MELTING_POINT_VIEW = {
    "Record": {
        "RecordNumber": 2244,
        "RecordTitle": "Aspirin",
        "Section": [
            {
                "TOCHeading": "Experimental Properties",
                "Section": [
                    {
                        "TOCHeading": "Melting Point",
                        "Information": [
                            {
                                "ReferenceNumber": 3,
                                "Value": {"Number": [135], "Unit": "°C"},
                            },
                            {"ReferenceNumber": 3, "Value": {"Number": [135.5, 136.0]}},
                        ],
                    }
                ],
            }
        ],
        "Reference": [{"ReferenceNumber": 3, "SourceName": "HSDB"}],
    }
}

VIEW_INDEX = {
    "Record": {
        "RecordNumber": 2244,
        "Section": [
            {
                "TOCHeading": "Names and Identifiers",
                "Section": [{"TOCHeading": "Computed Descriptors"}],
            },
            {
                "TOCHeading": "Safety and Hazards",
                "Section": [
                    {
                        "TOCHeading": "Hazards Identification",
                        "Section": [{"TOCHeading": "GHS Classification"}],
                    }
                ],
            },
        ],
    }
}

NOT_FOUND = fault("PUGREST.NotFound", "No data found")


def mock_pubchem(router: respx.Router) -> respx.Router:
    """Register the knowledge-tool endpoints for CID 2244 (aspirin)."""
    rx = re.escape
    router.post(
        url__regex=rx(REST) + r"/compound/fastsimilarity_2d/smiles/property/.*"
    ).respond(200, json=table(SALICYLIC, BROKEN, ASPIRIN, METHYL_SALICYLATE))
    router.post(
        url__regex=rx(REST) + r"/compound/fastsubstructure/smiles/property/.*",
        name="substructure",
    ).respond(200, json=table(METHYL_SALICYLATE, ASPIRIN))
    router.post(
        url__regex=rx(REST) + r"/compound/fastidentity/smiles/property/.*"
    ).respond(200, json=table(ASPIRIN))
    router.get(rx_url(REST + "/compound/cid/2244/synonyms/JSON")).respond(
        200, json=SYNONYMS
    )
    router.get(
        rx_url(REST + "/compound/cid/2244/xrefs/RegistryID,PubMedID,PatentID/JSON")
    ).respond(200, json=XREFS)
    router.get(rx_url(REST + "/compound/cid/2244/description/JSON")).respond(
        200, json=DESCRIPTION
    )
    router.get(rx_url(REST + "/compound/cid/2244/assaysummary/JSON")).respond(
        200, json=ASSAY_SUMMARY
    )
    for cid in (5, 6):  # untested / unlinked compounds
        for op in ("assaysummary", "xrefs/RegistryID,PubMedID,PatentID", "synonyms"):
            router.get(rx_url(f"{REST}/compound/cid/{cid}/{op}/JSON")).respond(
                404, json=NOT_FOUND
            )
    router.get(
        rx_url(VIEW + "/data/compound/2244/JSON"),
        params={"heading": "GHS Classification"},
    ).respond(200, json=GHS_VIEW)
    router.get(
        rx_url(VIEW + "/data/compound/2244/JSON"),
        params={"heading": "Experimental Properties"},
    ).respond(200, json=MELTING_POINT_VIEW)
    router.get(rx_url(VIEW + "/data/compound/2244/JSON")).respond(
        404, json=fault("PUGVIEW.NotFound", "No data found")
    )
    router.get(rx_url(VIEW + "/index/compound/2244/JSON")).respond(200, json=VIEW_INDEX)
    router.get(url__regex=rx(VIEW) + r"/(data|index)/compound/7/JSON.*").respond(
        404, json=fault("PUGVIEW.NotFound", "Invalid record number")
    )
    return router


def rx_url(url: str):
    """Exact URL pattern ignoring the query string."""
    return httpx.URL(url)
