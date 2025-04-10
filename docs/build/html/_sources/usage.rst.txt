How to Use BioMaster
====================

Using from Terminal
-------------------

1. Set API key and base URL
2. Move script into root folder
3. Download data to `./data/`
4. Set `id`
5. Run:

.. code-block:: bash

   python run.py

Results are saved in `./output/{id}/`.

Using the UI
------------

1. Run the UI:

.. code-block:: bash

   python runv.py

2. Visit http://127.0.0.1:7860/

3. Set parameters, define task, generate and execute plan interactively.

Output Format
-------------

- PLAN file: `output/{id}_PLAN.json`
- Scripts: `output/{id}_Step_{n}.sh`
- Logs: `output/{id}_DEBUG_Output_{n}.json`
