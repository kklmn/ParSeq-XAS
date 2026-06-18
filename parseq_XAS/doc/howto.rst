.. _howto:

ParSeq-XAS How-tos
------------------

Pipeline launch and command line options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ParSeq-XAS pipeline is started by running::

    python XAS_start.py

.. hint::

   Use the ``-h`` option to display the available startup parameters. Two
   particularly useful options are ``-p {filename}`` to load an existing
   project file and ``-v 100`` to increase verbosity for troubleshooting.

Data loading from the file tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. imagezoom:: _images/filemenu.png
   :align: right
   :loc: upper-right-corner
   :alt: &ensp;A popup menu over a column file in the file tree.

.. |icoLast| image:: _images/last.png
   :width: 12

Please see |formats|.

Use the |icoLast| button at the top of the file tree to return to the most
recently accessed file location and reuse the associated data format definition.

The content of a column file can be displayed directly within the ParSeq
application using the Metadata splitter. This option is available from the
context menu (right-click) in the file tree.

For beamlines that produce data files with variable formats depending on the
instruments used, ParSeq supports automation of format definitions -- provided
that the files include a header line describing the columns. The ParSeq-XAS
pipeline implements ``auto_format()`` methods in several data nodes (see the
``XAS_nodes`` module), which can be customized for specific requirements.
Auto-formatting can be triggered from the context menu.

Once a format is defined for the selected column file or HDF5 entry, ParSeq
attempts to import the data. If the import is successful, the selection is
highlighted in green; otherwise, it appears in red.

Data loading can then be initiated either from the context menu or by dragging
and dropping items from the file tree onto the data tree.

Data loading from an external file browser
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data files can be dragged and dropped from an external file browser onto the
data tree. If data sources are defined in the data format widget, they will be
applied to the dropped files. Otherwise, the ``auto_format()`` procedure
(see above) will be attempted.

Data range
~~~~~~~~~~

The data format widget includes a "conversion" tab. Refer to its tooltip panel
for detailed usage information. In particular, the function ``lim(Emin, Emax)``
can be used in the topmost input field to restrict the energy range to the
desired interval.

Data deglitching
~~~~~~~~~~~~~~~~

.. imagezoom:: _images/XAS-glitch.gif
   :align: right
   :alt: &ensp;A demonstration of glitch removal by scaling.

Please see |corrections|.

Although data corrections can be applied in any transformation node, removal of
monochromator glitches is most straightforward in the :math:`\chi(k)` node.

**Note**: When performing deglitching in the µd(E) node, ensure that pre-edge
subtraction and edge normalization are disabled.

.. raw:: html

   <div class="clearer"> </div>

Data combinations
~~~~~~~~~~~~~~~~~

The following data combinations are available: average, sum, RMS deviation, and
only for 1D: classical PCA, cumulative PCA, target transformation and
MCR-ALS. Please examine |combinations|.

Project files and data saving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParSeq-XAS includes several example project files located in the ``saved``
folder. Use the slider in the preview panel to browse through images of the
pipeline nodes:

.. imagezoom:: _images/load-proj.gif
   :align: center
   :alt: &ensp;Preview in a ParSeq project file.

.. note::

   Project files and their associated data files can usually be moved to a new
   location without losing the references stored in the project file. This may
   not work when files are located on network storage.

.. note::

   When saving a project file, pay attention to the current data selection.
   Only the selected data items will be exported.

Make publication plots
~~~~~~~~~~~~~~~~~~~~~~

.. imagezoom:: _images/save-proj.png
   :align: right
   :alt: &ensp;Saving a ParSeq project file.

1. ParSeq plot windows (based on silx plots) provide a Save button that can
   export the current plot view to a graphics format.

2. The Save Project dialog in ParSeq includes an option to generate a plotting
   script together with the exported data. These scripts contain commented
   sections for adjusting energy ranges, colors, and other plot settings to
   facilitate further customization.

3. The launch script ``XAS_start.py`` can be run with the ``-p`` option to load
   a project file and with the ``-nG`` option to execute the pipeline without
   a GUI. In this headless mode, the script includes a dedicated section that
   performs plotting using ``matplotlib``. This section can be freely modified
   as needed.
