import abc
try:
    # Python 3.x
    import collections.abc as collections
except ImportError:
    # Python 2.x
    import collections

from ..plugins.PyScript import DataObject, CloneHelper, DataCollection, AttributeDataObject

# Implementation of the DataCollection.attributes field.
def _DataCollection_attributes(self):
    """
    A dictionary of key-value pairs that represent global tokens of information which are not associated with
    any specific data object in the data collection. 
    
    An *attribute* is a value of type ``int``, ``float``, or ``str`` with a unique identifier name such 
    as ``"Timestep"`` or ``"ConstructSurfaceMesh.surface_area"``. The attribute name serves as keys for the :py:attr:`!attributes` dictionary of the data collection.
    Attributes are dynamically generated by modifiers in a data pipeline or by a data source as explained in the following.
    
    **Attributes loaded from input files**

    The ``Timestep`` attribute is loaded from LAMMPS dump files and other simulation file formats
    that store the simulation timestep. Such input attributes can be retrieved from 
    the :py:attr:`!.attributes` dictionary of a pipeline's :py:attr:`~ovito.pipeline.FileSource`::

        >>> pipeline = import_file('snapshot_140000.dump')
        >>> pipeline.source.attributes['Timestep']
        140000
        
    Other attributes read from an input file are, for example, the key-value pairs found in the header line of *extended XYZ* files.
    
    **Dynamically computed attributes**
    
    Analysis modifiers like the :py:class:`~ovito.modifiers.CommonNeighborAnalysisModifier` or
    the :py:class:`~ovito.modifiers.ClusterAnalysisModifier` output scalar computation results
    as attributes. The reference documentation of each modifier type lists the attributes it produces.
    
    For example, the number of clusters identified by the :py:class:`~ovito.modifiers.ClusterAnalysisModifier`
    can be queried as follows::
    
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff = 3.1))
        data = pipeline.compute()
        nclusters = data.attributes["ClusterAnalysis.cluster_count"]
        
    **Exporting attributes to a text file**
    
    The :py:func:`ovito.io.export_file` function supports writing attribute values to a text
    file, possibly as functions of time::
    
        export_file(pipeline, "data.txt", "txt", 
            columns = ["Timestep", "ClusterAnalysis.cluster_count"], 
            multiple_frames = True)
            
    **User-defined attributes**
    
    The :py:class:`~ovito.modifiers.PythonScriptModifier` allows you to generate your own
    attributes that are dynamically computed (typically on the basis of some other input information):
    
    .. literalinclude:: ../example_snippets/python_modifier_generate_attribute.py
        :lines: 6-

    The :py:class:`~ovito.modifiers.CommonNeighborAnalysisModifier` used in the example above generates
    the attribute ``CommonNeighborAnalysis.counts.FCC`` to report the number of atoms that 
    form an FCC lattice. To compute the fraction of FCC atoms from that, we need to divide by the total number of 
    atoms in the system. To this end, we insert a :py:class:`~ovito.modifiers.PythonScriptModifier` 
    into the pipeline behind the :py:class:`~ovito.modifiers.CommonNeighborAnalysisModifier`.
    Our custom modifier function generates a new attribute named ``fcc_fraction``. Finally, 
    the value of the user-defined attribute can be queried from the pipeline or exported to 
    a text file using the :py:func:`~ovito.io.export_file` function as described above.
    """
        
    # Helper class used to implement the DataCollection.attributes field.
    class _AttributesView(collections.MutableMapping):
        
        def __init__(self, data_collection):
            """ Constructor that stores away a back-pointer to the owning DataCollection instance. """
            self._collection = data_collection
            
        def __len__(self):
            count = 0
            for obj in self._collection.objects:
                if isinstance(obj, AttributeDataObject):
                    count += 1
            return count
        
        def __getitem__(self, key):
            for obj in self._collection.objects:
                if isinstance(obj, AttributeDataObject) and obj.id == key:
                    return obj.value
            raise KeyError("Attribute '%s' does not exist in data collection." % key)

        def __setitem__(self, key, value):
            for obj in self._collection.objects:
                if isinstance(obj, AttributeDataObject) and obj.id == key:
                    if not value is None:
                        obj.value = value
                    else:
                        del self._collection.objects[obj]
                    return
            if not value is None:
                attr =  AttributeDataObject(id = key, value = value)
                self._collection.objects.append(attr)

        def __delitem__(self, key):
            """ Removes a global attribute from the data collection. """
            for obj in self._collection.objects:
                if isinstance(obj, AttributeDataObject) and obj.id == key:
                    del self._collection.objects[obj]
                    return
            raise KeyError("Attribute '%s' does not exist in data collection." % key)

        def __iter__(self):
            """ Returns an iterator over the names of all global attributes. """
            for obj in self._collection.objects:
                if isinstance(obj, AttributeDataObject):
                    yield obj.id

        def __repr__(self):
            return repr(dict(self))
    
    return _AttributesView(self)
DataCollection.attributes = property(_DataCollection_attributes)

# Implementation of the DataCollection.copy_if_needed() method.
def _DataCollection_copy_if_needed(self, obj, deepcopy=False):
    """
    Makes a copy of a data object from this data collection if the object is not exclusively 
    owned by the data collection but shared with other collections. After the method returns,
    the data object is exclusively owned by the collection and it becomes safe to modify the object without
    causing unwanted side effects.

    Typically, this method is used in custom modifier functions (see :py:class:`~ovito.modifiers.PythonScriptModifier`) that
    participate in OVITO's data pipeline system. A modifier function receives an input collection of
    data objects from the system. However, modifying these input
    objects in place is not allowed, because they are owned by the pipeline and modifying them would 
    lead do unexpected side effects.
    This is where this method comes into play: It makes a copy of a given data object and replaces
    the original in the data collection with the copy. The caller can now safely modify this copy in place,
    because no other data collection can possibly be referring to it.

    The :py:meth:`!copy_if_needed` method first checks if *obj*, which must be a data object from this data collection, is
    shared with some other data collection. If yes, it creates an exact copy of *obj* and replaces the original
    in this data collection with the copy. Otherwise it leaves the object as is, because it is already exclusively owned
    by this data collection. 

    :param DataObject obj: The object from this data collection to be copied if needed.
    :return: An exact copy of *obj* if it was shared with some other data collection. Otherwise the original object is returned.
    """
    assert(isinstance(obj, DataObject))
    # The object to be modified must be in this data collection.
    if obj not in self.objects:
        raise ValueError("DataCollection.copy_if_needed() must be called with an object that is part of the data collection.")
    if obj.num_strong_references <= 1:
        return obj
    idx = self.objects.index(obj)
    clone = CloneHelper().clone(obj, deepcopy)
    self.objects[idx] = clone
    assert(clone.num_strong_references == 1)
    return clone
DataCollection.copy_if_needed = _DataCollection_copy_if_needed

# Implementation of the DataCollection.find() method.
def _DataCollection_find(self, object_type):
    """
    Looks up the first data object from this collection of the given class type.

    :param object_type: The :py:class:`DataObject` subclass that should be looked up.
    :return: The first instance of the given class or its subclasses from the :py:attr:`.objects` list; or ``None`` if there is no instance.

    Method implementation::

        def find(self, object_type):
            for o in self.objects:
                if isinstance(o, object_type): return o
            return None
    """
    if not issubclass(object_type, DataObject):
        raise ValueError("Not a subclass of ovito.data.DataObject: {}".format(object_type))
    for obj in self.objects:
        if isinstance(obj, object_type):
            return obj
    return None
DataCollection.find = _DataCollection_find

# Implementation of the DataCollection.expect() method.
def _DataCollection_expect(self, object_type):
    """
    Looks up the first data object in this collection of the given class type.
    Raises a ``KeyError`` if there is no instance matching the type. Use :py:meth:`.find` instead
    to test if the data collection contains the given type of data object.

    :param object_type: The :py:class:`DataObject` subclass specifying the type of object to find.
    :return: The first instance of the given class or its subclasses from the :py:attr:`.objects` list.
    """
    o = self.find(object_type)
    if o is None: raise KeyError('Data collection does not contain a {}'.format(object_type))
    return o
DataCollection.expect = _DataCollection_expect

# Implementation of the DataCollection.find_all() method.
def _DataCollection_find_all(self, object_type):
    """
    Looks up all data objects from this collection of the given class type.

    :param object_type: The :py:class:`DataObject` subclass that should be looked up.
    :return: A Python list containing all instances of the given class or its subclasses from the :py:attr:`.objects` list.

    Method implementation::

        def find_all(self, object_type):
            return [o for o in self.objects if isinstance(o, object_type)]
    """
    if not issubclass(object_type, DataObject):
        raise ValueError("Not a subclass of ovito.data.DataObject: {}".format(object_type))
    return [obj for obj in self.objects if isinstance(obj, object_type)]
DataCollection.find_all = _DataCollection_find_all

# This method has been deprecated and is here only for backward compatibility with OVITO 2.9.0.
# Use DataCollection.objects.append() instead.
def _DataCollection_add(self, obj):
    if not obj in self.objects:
        self.objects.append(obj)
DataCollection.add = _DataCollection_add

# This method has been deprecated and is here only for backward compatibility with OVITO 2.9.0.
# Use del DataCollection.objects[] instead.
def _DataCollection_remove(self, obj):
    index = self.objects.index(obj)
    if index >= 0:
        del self.objects[index]
DataCollection.remove = _DataCollection_remove

# This method has been deprecated and is here only for backward compatibility with OVITO 2.9.0.
def _DataCollection_replace(self, oldobj, newobj):
    index = self.objects.index(oldobj)
    if index >= 0:
        self.objects[index] = newobj
DataCollection.replace = _DataCollection_replace
