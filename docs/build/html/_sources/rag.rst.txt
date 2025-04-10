Knowledge Management with RAG
=============================

PLAN RAG
--------

To add a new analysis workflow:

1. Collect reference workflow
2. Write it using standard format
3. Add entry to `doc/Plan_Knowledge.json`

EXECUTE RAG
-----------

1. Document tool/script/function usage
2. Add to `doc/Task_Knowledge.json`

Use `source` metadata for retrieval.

Maintain quality by:
- Being specific
- Avoiding redundant entries
- Deleting `chroma_db` after changes
