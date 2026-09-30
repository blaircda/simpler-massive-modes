---
name: hepth-arxiv-inspire-literature-exploration
description: >-
  Use this skill during literature search or topic research in theoretical high-energy
  physics.
---
# HEP Literature Traversal and Citation Searching
When conducting research on a topic and encountering a promising paper on arXiv:
## 1. Retrieve the Seed Paper on INSPIRE-HEP
Extract the arXiv identifier (`YYMM.NNNNN` or old format `hep-th/XXXXXXX`) and fetch its record:
- **API Endpoint**: `https://inspirehep.net/api/literature?q=arxiv:<ARXIV_ID>`
- **Capture**:
  - Record ID: `recid = hits.hits[0].id`
  - Core metadata: authors, publication date, abstract, keywords
  - Total citation count (`citation_count`)
## 2. Follow-ups
Query papers that cite the seed paper to find how the idea evolved:
- **Most Recent Follow-ups**:
  `https://inspirehep.net/api/literature?q=refersto%3Arecid%3A<RECID>&sort=mostrecent&size=50&page=1`
  `https://inspirehep.net/api/literature?q=refersto%3Arecid%3A<RECID>&sort=mostrecent&size=50&page=2`
  and so on
- **Most Influential Extensions**:
  `https://inspirehep.net/api/literature?q=refersto%3Arecid%3A<RECID>&sort=mostcited&size=25`
## 4. Relevance Filtering Heuristics
You can adjust the size parameter in the above query for a more focused selection, however in some cases what we want to find might be hidden in an old obscure paper with few citations. This will depnd on the topic.
In any case, do not read all papers in full indiscrimately. 
Filter using these heuristics:
1. **Topic Alignment**: Scan titles and abstracts for specific technical keywords relevant to the topic.
2. **Author Continuity**: Look for follow-ups by the same authors or their known collaborators (often indicates the next paper in a research program).
## 5. Subagent Research Synthesis
When reporting back to the parent agent or user, make sure to provide the explicit inspirehep.net URL for the papers you have identified as being of interest.
