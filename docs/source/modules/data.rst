Data
====
This is sub-packege for testing. It has ``ligands``, ``receptors``, ``boxes`` and ``config_yaml``.
There are two sytems: r_x0161 and r_6lu7.  The corresponded atributes could be acces
like this: 

.. ipython:: python

    from moldrug.data import ligands
    print(ligands.r_x0161)
    print(ligands.r_6lu7)

.. ipython:: python

    from moldrug.data import boxes
    print(boxes.r_x0161)
    print(boxes.r_6lu7)