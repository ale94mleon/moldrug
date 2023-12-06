Data
====
This is a sub-package for testing. It has the function ``get_data`` to retrieve the data.

.. ipython:: python

    from moldrug.data import get_data
    import json
    print(json.dumps(get_data('x0161'), indent = 3))
    print(json.dumps(get_data('6lu7'), indent = 3))