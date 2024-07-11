=======================================
Additional Information
=======================================

The following sections contains information on some suplimental scripts provided with HypOptLib, as well as general notes on
the implementation.

Analysis Scripts
=======================================

Several Python scripts are provided for basic data analysis. They are intended mostly as a guide for how to work with HypOptLib
data. the intent is mostly to provide a starting point for developers own analysis.

utilities.py
---------------------------------------

The utilities file contains some helpfer functions, but the main object is the xdmf helper class.
This class provides a basic xmf writer to interface ParaView and other viewers with the HypOptLib
HDF5 output files.

Methods
_______________________________________
startGenericFile()
    Writes the header for an arbitrary xmf file. Must be closed with closeGenericFile().

closeGenericFile()
    Writes the footer to a generic xmf file.

insertRectilinearMesh(dimensions)
    Inserts a rectilinear mesh of the provided dimensions.

insertAttribute(fileName, attributeName, attributeLocation)
    Inserts the desired attribute into the xmf file.

startTimesteppedFile()
    Writes the header for a timestepped xmf file. Must be closed with closeTimesteppedFile().

closeTimesteppedFile()
    Writes the footer to a timestepped xmf file.

insertTimestep(time, hdf5FileName, iterationName)
    Inserts a timestep into an opened timestepped xmf file.

animate.py
---------------------------------------

The animate script generates an xdmf file with all timestep information. This allows ParaView
to properly interpret the HDF5 output file as a timestepped file, and in turn create animations.
It supports an arbitrary number of restart files.


finalPositionToXdmf.py
---------------------------------------

This script creates a simple xdmf file pointing to the final positions of the provided file. This is mostly
designed to show how to access iterations in ParView.

generateMeanMap.py
---------------------------------------

This script calculates the mean of the position from the provided files, and writes it to a new xdmf and associated
HDF5 file. This is a useful script to show how a system and structure changes behaves at different temperatures.

A Note on parallelization
=======================================
Hyperoptimization itself is not an algorithm that is inherently parallelisable. There are a few steps that must be sequential,
or that must break parallelism. For example, there are a number of sections that rely on left-shifting by one element, and reducing
the size of that vector. This is mostly related to Nose Hoover calculations. This means that while the rest of the algorithm can be done
in parallel, in general it should be considered a sequential algorithm.

However, this limitation does not prevent using parallel techniques when solving the objective function, sensitivities, filtering, etc.
In fact, HypOptLb has been designed to take advantage of parallel computing to facilitate these operations. For example, the finite
element analysis done in the Linear Elasticity class from TopOpt uses a parallelized KSP solver from PETSc. So, since HypOptLib has been
designed in order to support mpi parallelism, users can take full advantage of the parallel KSP solvers to greatly reduce the time to solve
problems. This is why we refer to HypOptLib as semi-parallel: there are sequential parts of the Hyperoptimization algorithm, but the rest is
implemented to support parallel computing in order to allow for faster computation of the finite element analysis.

File Formats
=======================================
HypOptLib uses HDF5 for the simulation save files, and for all telemetry and restart information. This could be modified if desired to use other
file formats supported by PETSc, such as raw binary. All the file access is handled by the FileManager class, so any modifications to file type
or format can be done there.

The data structure is designed as follows:

.. code-block:: text

    file.h5
    |___Dataset
    |       Compliance
    |       Energy Error
    |       FEA Solver Itr Hamiltonian
    |       FEA Solver Itr Sensitivity
    |       Final Position
    |       Final State Field
    |       Final Velocity
    |       Final even NH Pos
    |       Final even NH Vel
    |       Final odd NH Pos
    |       Final odd NH Vel
    |       Hamiltonian
    |       Iteration Compute Time
    |       Lambda
    |
    |_______State
    |       |   iteration0
    |       |   iteration1
    |       |   ...
    |
    |       Temperature
    |       Timestep
    |       Volume Fraction
    |
    |___Setting
    |       Init Conditions File
    |       Random Init Conditions
    |       T
    |       ch_order
    |       dt
    |       nelx
    |       nely
    |       nelz
    |       penal
    |       rmin
    |       saveFreq
    |       stepN
    |       volfrac
    |
    |_______boundaryConditions
    |       |   boundaryCondition_0
    |       |   boundaryCondition_1
    |       |   ...
    |
    |_______domain


Visualization and ParaView
=======================================
ParaView is the recommended tool for visualising HypOptLib outputs. The analysis scripts mentioned above can be used to generate xmf files that
ParaView can read both for animations and static models. The following section outlines some general tips for those unfamiliar with ParaView.

Viewing Volumes
---------------------------------------
After generating the xmf file, open the file in ParaView. You will be prompeted to select a reader with which to view the file. Usually the default
one is fine.

.. figure:: ../figures/paraview_open_file_prompt.png
    :width: 500
    :align: center

The file with then be loaded in the Pipeline Browser. If there are multiple datasets present in the xmf file, select the one you'd like to view.

.. figure:: ../figures/paraview_press_apply.png
    :width: 500
    :align: center

Press **Apply** to load it in.

.. figure:: ../figures/paraview_select_volume.png
    :width: 800

Select the box that says **Outline**, and choose **Volume** from the dropdown.

.. figure:: ../figures/paraview_volume_selected.png
    :width: 800

The volume will then be displayed in the Layout window. If this is a timestepped file, then it will show the first timestep in the file.

.. figure:: ../figures/paraview_volume_open.png
    :width: 800

Animating Timestepped Files
---------------------------------------

Timestepped files can be animated through the **Time Inspector**. Under View, make sure this option is selected:

.. figure:: ../figures/paraview_time_inspector.png
    :width: 500
    :align: center

Then select the **Time Inspector** tab on the right of the Layout window:

.. figure:: ../figures/paraview_time_opened.png
    :width: 800

You can now press the play button to view the animation, or skip forward and backwards with the arrows.

For more advanced control, select the **Animation View** from the View dropdown.

.. figure:: ../figures/paraview_animation_view.png
    :width: 800

Different modes can be selected from the **Mode** drop down:

 * **Sequence** allows you to choose exactly how many frames to use in the animation.

 * **Real Time** uses all frames but selects the desired length of the animation

 * **Snap To Timesteps** uses every timestep as a frame.

To save the animatino, go to File and select **Save Animation**.

.. figure:: ../figures/paraview_save_animation.png
    :width: 800

After choosing a name and path to save the animation, you will be given options for the animation. The resolution can
be update here, along with frame rate an window to save. 

.. figure:: ../figures/paraview_animation_options.png
    :width: 500
    :align: center

This will render the animation and save it to the selected location.

.. figure:: ../figures/animation_example.gif
    :width: 700
    :align: center

Filters
---------------------------------------

There are many filters ParaView provides to analyze data, but as a basic introduction here is an introcution to a filter 
is especially useful for topology problems with HypOptLib.

The threshold filter allows you to filter out cells within an upper and lower bound. The main use for this
filter is to understand the internal structure of a simulation. For example, we will walk through looking at
thresholding a beam at temperature 0 early in the simulation. 

To select a filter, go to the Filters dropdown and either search for the desired filter, or browse the available categories.

.. figure:: ../figures/paraview_filters.png
    :width: 500
    :align: center

After selecting the Threshold filter, a new Pipeline will appear in the Pipeline Browser.

.. figure:: ../figures/paraview_threshold_browser.png
    :width: 500
    :align: center

Here you can select upper and lower threshold values, then press **Apply**. If the threshold is set to 0, then all cells with
a value of 0 will show as part of the surface, and the structure will not be visible.

.. figure:: ../figures/paraview_threshold_0.png
    :width: 800

However, as the lower threshold is increased, more of the internal structure is revealed. At threshold 1e-09, we can see the shape
of the beam start to be revealed. However, it is mostly covered in low density points.

.. figure:: ../figures/paraview_threshold_small.png
    :width: 800

At threshold 0.2, the lower density points are pulled back to reveal more structure.

.. figure:: ../figures/paraview_threshold_med.png
    :width: 800

And finally at threshold 0.6, msot of the internal structure is revealed, as the beam starts to take shape.

.. figure:: ../figures/paraview_threshold_large.png
    :width: 800


The other case in which the threshold is useful is to render high-resolution images, as the **Volume** view is not always as good
looking as the **Surface** view. For example, the following image shows the same simulation as the above example, but much further
in the simulation. This uses a threshold of 1e-09.

.. figure:: ../figures/paraview_threshold_small_end.png
    :width: 800

Notice how much clearer the above image is than the volume render below:

.. figure:: ../figures/paraview_threshold_volume_compare.png
    :width: 800

While the volume render appears totally converged, the thresholded image clearly shows lots of cells with low density, indicating the
simulation still needs more time to fully converge.

There are many more filters provided by ParaView which may be helpful. For more information, please see the ParaView Documentation
`here <https://docs.paraview.org/en/latest/UsersGuide/filteringData.html>`_.
